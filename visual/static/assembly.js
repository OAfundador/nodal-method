'use strict';
(() => {
  const $ = id => document.getElementById(id);
  const clone = value => JSON.parse(JSON.stringify(value));
  const fmt = (value, digits=3) => Number(value).toLocaleString('pt-BR', {minimumFractionDigits:digits,maximumFractionDigits:digits});
  const empty = () => ({version:1,kind:'assembly',name:'Associação de componentes',coupling:'independent',components:{},instances:[],view:{axis:'y',offset_m:0,region:'',unit:'C'}});
  const svgNode = (tag, attrs={}, text=null) => {const el=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [key,value] of Object.entries(attrs))el.setAttribute(key,value);if(text!==null)el.textContent=text;return el;};
  const option = (value, label) => {const el=document.createElement('option');el.value=value;el.textContent=label;return el;};
  class AssemblyApp {
    constructor(){this.data=empty();this.result=null;this.history=[];this.selected=0;this.busy=false;this.offline=false;this.editors=new Map();this.bind();}
    async post(route,payload){const response=await fetch('/api/assembly/'+route,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(payload)});const data=await response.json();if(!response.ok)throw new Error(data.error||'Falha na associação.');return data;}
    toast(text){$('toast').textContent=text;$('toast').hidden=false;clearTimeout(this.toastTimer);this.toastTimer=setTimeout(()=>$('toast').hidden=true,7000);}
    change(action){if(this.busy)return;const before=clone(this.data);try{action();this.history.push(before);if(this.history.length>35)this.history.shift();this.result=null;this.render();}catch(error){this.data=before;this.render();this.toast(error.message);}}
    load(data){this.change(()=>{this.data=clone(data);this.data.view={...empty().view,...this.data.view};this.selected=0;});}
    key(prefix,used){let i=1;while(used.includes(prefix+i))i++;return prefix+i;}
    download(name,content,type){const url=URL.createObjectURL(new Blob([content],{type})),link=document.createElement('a');link.href=url;link.download=name;link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
    async init(){
      const saved=$('saved-assembly');
      if(saved){try{const data=JSON.parse(saved.textContent);this.data=data.assembly;this.result=data.result;this.offline=true;document.body.classList.add('offline');}catch(error){this.toast('O relatório salvo está incompleto.');}}
      this.render();
      if(!this.offline&&window.opener)window.opener.postMessage({type:'assembly-ready'},location.origin);
    }
    bind(){
      $('open-assembly').addEventListener('click',()=>$('assembly-file').click());
      $('add-component').addEventListener('click',()=>$('component-file').click());
      $('assembly-file').addEventListener('change',async event=>{try{const file=event.target.files[0];if(!file)return;if(file.size>1048576)throw new Error('Use um arquivo de até 1 MB.');const text=await file.text(),data=await this.post(file.name.toLowerCase().endsWith('.txt')?'import':'validate',file.name.toLowerCase().endsWith('.txt')?{text}:JSON.parse(text));this.load(data);this.toast('Associação aberta.');}catch(error){this.toast(error.message);}finally{event.target.value='';}});
      $('component-file').addEventListener('change',async event=>{try{for(const file of event.target.files){if(file.size>1048576)throw new Error('Use componentes de até 1 MB.');await this.addComponent(JSON.parse(await file.text()));}}catch(error){this.toast(error.message);}finally{event.target.value='';}});
      $('load-example').addEventListener('click',async()=>{try{const response=await fetch('/api/assembly/example');if(!response.ok)throw new Error('Exemplo indisponível.');this.load(await response.json());}catch(error){this.toast(error.message);}});
      $('assembly-name').addEventListener('change',()=>this.change(()=>this.data.name=$('assembly-name').value||'Associação'));
      $('add-instances').addEventListener('click',()=>this.addInstances());
      $('undo-assembly').addEventListener('click',()=>{if(this.busy||!this.history.length)return;this.data=this.history.pop();this.result=null;this.render();});
      $('solve-assembly').addEventListener('click',()=>this.solve());
      $('save-json').addEventListener('click',async()=>{try{const data=this.offline?this.data:await this.post('validate',this.data);this.download('associacao.json',JSON.stringify(data,null,2),'application/json');}catch(error){this.toast(error.message);}});
      $('save-txt').addEventListener('click',async()=>{try{const data=await this.post('repl',this.data);this.download('associacao.txt',data.text,'text/plain;charset=utf-8');}catch(error){this.toast(error.message);}});
      $('assembly-csv').addEventListener('click',()=>this.exportCSV());$('assembly-png').addEventListener('click',()=>this.exportPNG());$('assembly-report').addEventListener('click',()=>this.exportReport());
      $('selected-instance').addEventListener('change',()=>{this.selected=Number($('selected-instance').value);this.renderResults();});
      for(const id of ['view-axis','view-region','view-offset','view-unit'])$(id).addEventListener('change',()=>{const offset=Number($('view-offset').value)/1000;if(!Number.isFinite(offset)){this.toast('Informe uma origem numérica.');return;}this.data.view={axis:$('view-axis').value,region:$('view-region').value,offset_m:offset,unit:$('view-unit').value};this.renderResults();});
      $('assembly-map').addEventListener('click',event=>{const cell=event.target.closest('[data-instance]');if(cell){this.selected=Number(cell.dataset.instance);this.renderResults();}});
      $('assembly-map').addEventListener('pointermove',event=>{const cell=event.target.closest('[data-instance]');if(cell)$('map-readout').textContent=cell.getAttribute('aria-label');});
      window.addEventListener('message',async event=>{
        if(this.offline||event.origin!==location.origin)return;
        try{
          if(event.data?.type==='component-draft'&&event.source===window.opener)await this.addComponent(event.data.model);
          if(event.data?.type==='editor-ready'&&this.editors.has(event.source)){const key=this.editors.get(event.source);event.source.postMessage({type:'edit-model',model:this.data.components[key],componentKey:key},location.origin);}
          if(event.data?.type==='component-updated'&&this.editors.get(event.source)===event.data.componentKey)await this.addComponent(event.data.model,event.data.componentKey);
        }catch(error){this.toast(error.message);}
      });
    }
    async addComponent(model,key=null){
      if(model?.kind==='assembly')throw new Error('Este arquivo contém uma associação. Use Abrir associação.');
      key=key||this.key('componente',Object.keys(this.data.components));
      const proposal=clone(this.data);proposal.components[key]=clone(model);await this.post('validate',proposal);
      this.change(()=>this.data.components[key]=clone(model));$('repeat-component').value=key;this.toast('Componente disponível para reutilização.');
    }
    addInstances(){this.change(()=>{
      const component=$('repeat-component').value;if(!this.data.components[component])throw new Error('Importe um componente primeiro.');
      const count=Number($('repeat-count').value),text=$('repeat-factors').value.trim(),factors=text?text.split(/[\s;,]+/).map(Number):[];
      if(!Number.isInteger(count)||count<1||this.data.instances.length+count>100)throw new Error('Use até 100 instâncias por associação.');
      if(factors.length&&factors.length!==count)throw new Error('Informe exatamente um fator por instância, ou deixe a lista vazia.');
      if(factors.some(v=>!Number.isFinite(v)||v<0))throw new Error('Os fatores devem ser números maiores ou iguais a zero.');
      const prefix=$('repeat-prefix').value.trim()||'Instância',start=this.data.instances.length;
      for(let i=0;i<count;i++){const id=this.key('inst',this.data.instances.map(v=>v.id));this.data.instances.push({id,name:prefix+' '+(start+i+1),component,source_factor:factors.length?factors[i]:1});}
    });}
    async solve(){
      if(this.busy)return;this.busy=true;document.body.classList.add('busy');$('solve-assembly').disabled=true;$('assembly-status').textContent='Calculando as instâncias…';
      try{this.result=await this.post('solve',this.data);this.selected=Math.min(this.selected,this.result.instances.length-1);this.toast('Associação calculada. Clique nas colunas para inspecionar os componentes.');}
      catch(error){this.result=null;this.toast(error.message);}finally{this.busy=false;document.body.classList.remove('busy');$('solve-assembly').disabled=false;this.render();}
    }
    render(){
      document.title=(this.data.name||'Associação')+' · Método Nodal';$('assembly-name').value=this.data.name;
      const list=$('component-list');list.replaceChildren();const select=$('repeat-component'),previous=select.value;select.replaceChildren();
      for(const [key,model] of Object.entries(this.data.components)){
        select.append(option(key,model.name||key));const card=document.createElement('div');card.className='component-card';const name=document.createElement('strong');name.textContent=model.name||key;const info=document.createElement('p');info.textContent=model.regions.length+' regiões · '+model.mesh.nx+' × '+model.mesh.ny+' divisões · '+this.data.instances.filter(v=>v.component===key).length+' instâncias';
        const actions=document.createElement('div');actions.className='component-actions';
        for(const [label,action] of [['Editar no desenho',()=>{const editor=window.open('/','_blank');if(editor)this.editors.set(editor,key);else this.toast('Permita abrir a aba de edição.');}],['Duplicar',()=>this.addComponent({...clone(model),name:(model.name||key)+' · cópia'}).catch(e=>this.toast(e.message))],['Excluir',()=>this.change(()=>{if(this.data.instances.some(v=>v.component===key))throw new Error('Remova primeiro as instâncias que usam este componente.');delete this.data.components[key];})]]){const button=document.createElement('button');button.textContent=label;button.addEventListener('click',action);actions.append(button);}card.append(name,info,actions);list.append(card);
      }
      if(this.data.components[previous])select.value=previous;
      $('empty-library').hidden=!!Object.keys(this.data.components).length;
      const body=$('instance-list');body.replaceChildren();
      this.data.instances.forEach((instance,index)=>{
        const tr=document.createElement('tr');
        const name=document.createElement('input');name.value=instance.name;name.setAttribute('aria-label','Nome da instância '+(index+1));name.addEventListener('change',()=>this.change(()=>instance.name=name.value||instance.id));
        const component=document.createElement('select');component.setAttribute('aria-label','Componente da instância '+(index+1));for(const [key,model] of Object.entries(this.data.components))component.append(option(key,model.name||key));component.value=instance.component;component.addEventListener('change',()=>this.change(()=>instance.component=component.value));
        const factor=document.createElement('input');factor.type='number';factor.step='any';factor.min=0;factor.value=instance.source_factor;factor.setAttribute('aria-label','Fator de geração da instância '+(index+1));factor.addEventListener('change',()=>this.change(()=>{const value=Number(factor.value);if(factor.value===''||!Number.isFinite(value)||value<0)throw new Error('Informe um fator numérico maior ou igual a zero.');instance.source_factor=value;}));
        const remove=document.createElement('button');remove.textContent='Excluir';remove.setAttribute('aria-label','Excluir instância '+(index+1));remove.addEventListener('click',()=>this.change(()=>this.data.instances.splice(index,1)));
        for(const child of [name,component,factor,remove]){const td=document.createElement('td');td.append(child);tr.append(td);}body.append(tr);
      });
      $('instance-empty').hidden=!!this.data.instances.length;$('assembly-status').textContent=this.result?'Cálculo concluído · instâncias independentes':this.data.instances.length+' instâncias associadas.';$('undo-assembly').disabled=!this.history.length;
      this.renderResults();
    }
    unit(){return this.data.view.unit==='K'?'K':'°C';}
    temperature(value){return this.data.view.unit==='K'?value+273.15:value;}
    samples(instance,region=this.data.view.region){
      const axis=this.data.view.axis,groups=new Map();
      for(const cell of instance.result.cells){if(region&&cell.region!==region)continue;const lo=cell[axis+'0'],hi=cell[axis+'1'],key=lo+'|'+hi,old=groups.get(key);if(!old||cell.T_C>old.T_C)groups.set(key,{lo,hi,T_C:cell.T_C});}
      return Array.from(groups.values()).sort((a,b)=>a.lo-b.lo);
    }
    color(value,low,high){const t=high-low<1e-10?.5:Math.max(0,Math.min(1,(value-low)/(high-low))),stops=[[255,255,204],[254,217,118],[253,141,60],[227,26,28],[128,0,38]],p=t*4,i=Math.min(3,Math.floor(p)),a=p-i;return 'rgb('+stops[i].map((v,k)=>Math.round(v+(stops[i+1][k]-v)*a)).join(',')+')';}
    renderResults(){
      const view=this.data.view;$('view-axis').value=view.axis;$('view-offset').value=+(view.offset_m*1000).toPrecision(12);$('view-unit').value=view.unit;
      const regionSelect=$('view-region');regionSelect.replaceChildren(option('','Todas · máxima em cada faixa'));const regions=new Map();for(const model of Object.values(this.data.components))for(const r of model.regions)if(!regions.has(r.id))regions.set(r.id,r.name||r.id);for(const [key,name] of regions)regionSelect.append(option(key,name+' ['+key+']'));regionSelect.value=regions.has(view.region)?view.region:'';
      if(view.region&&!regions.has(view.region))view.region='';
      const ready=!!this.result;$('result-empty').hidden=ready;$('result-content').hidden=!ready;for(const id of ['assembly-csv','assembly-png','assembly-report'])$(id).disabled=!ready;
      $('result-title').textContent=ready?this.data.name:'3. Resultados do conjunto';if(!ready)return;
      $('result-count').textContent=this.result.count;$('result-maximum').textContent=fmt(this.temperature(this.result.maximum_C),4)+' '+this.unit();$('result-residual').textContent=this.result.max_residual_W.toExponential(2)+' W';$('result-cells').textContent=this.result.total_cells;
      const selected=$('selected-instance');selected.replaceChildren();this.result.instances.forEach((v,i)=>selected.append(option(i,v.name)));selected.value=this.selected;
      this.drawMap();this.drawProfile();
    }
    drawMap(){
      const svg=$('assembly-map');svg.replaceChildren();const instances=this.result.instances,series=instances.map(v=>this.samples(v)),all=series.flat();if(!all.length){svg.append(svgNode('text',{x:550,y:300,'text-anchor':'middle'},'Nenhuma célula corresponde à região selecionada.'));return;}
      const low=Math.min(...all.map(v=>v.T_C)),high=Math.max(...all.map(v=>v.T_C)),axis=this.data.view.axis,length=Math.max(...instances.map(v=>v.domain[axis==='y'?'height':'width'])),offset=this.data.view.offset_m,w=900/instances.length;
      svg.append(svgNode('text',{x:85,y:30,'font-size':16,fill:'#304f62'},(this.data.view.region?$('view-region').selectedOptions[0].textContent:'Todas as regiões')+' · máxima por faixa'),svgNode('text',{x:85,y:51,'font-size':10,fill:'#748592'},'Cada coluna é uma instância independente; coordenadas reais no eixo '+axis+'.'));
      series.forEach((samples,index)=>{for(const sample of samples){const label=instances[index].name+' · '+axis+' = '+fmt(((sample.lo+sample.hi)/2+offset)*1000)+' mm · '+fmt(this.temperature(sample.T_C),4)+' '+this.unit();const rect=svgNode('rect',{x:85+index*w,y:480-sample.hi/length*400,width:w+.15,height:(sample.hi-sample.lo)/length*400+.15,fill:this.color(sample.T_C,low,high),'data-instance':index,'aria-label':label});rect.append(svgNode('title',{},label));svg.append(rect);}
        if(index===this.selected)svg.append(svgNode('rect',{x:85+index*w,y:80,width:w,height:400,fill:'none',stroke:'#263f50','stroke-width':2,'pointer-events':'none'}));
        svg.append(svgNode('text',{x:85+(index+.5)*w,y:503,'text-anchor':'middle','font-size':instances.length>50?7:10,fill:'#617d8e'},index+1));});
      for(let i=0;i<=5;i++)svg.append(svgNode('text',{x:73,y:484-i*80,'text-anchor':'end','font-size':11,fill:'#718795'},fmt((length*i/5+offset)*1000,2)));
      svg.append(svgNode('text',{x:85,y:70,'font-size':10,fill:'#718795'},axis+' [mm]'),svgNode('text',{x:535,y:534,'text-anchor':'middle','font-size':12,fill:'#617d8e'},'Posição da instância na associação'));
      const defs=svgNode('defs'),gradient=svgNode('linearGradient',{id:'assembly-gradient'});for(let i=0;i<=16;i++)gradient.append(svgNode('stop',{offset:i/16,'stop-color':this.color(low+(high-low)*i/16,low,high)}));defs.append(gradient);svg.append(defs,svgNode('rect',{x:250,y:564,width:600,height:12,fill:'url(#assembly-gradient)'}),svgNode('text',{x:250,y:598,'font-size':12,fill:'#617d8e'},fmt(this.temperature(low),4)+' '+this.unit()),svgNode('text',{x:850,y:598,'text-anchor':'end','font-size':12,fill:'#617d8e'},fmt(this.temperature(high),4)+' '+this.unit()));
    }
    drawProfile(){
      const instance=this.result.instances[this.selected];if(!instance)return;$('selected-title').textContent=instance.name+' · '+instance.component+' · fator '+fmt(instance.source_factor,4);
      const svg=$('instance-profile');svg.replaceChildren();const series=instance.regions.map((r,i)=>({name:r.name,values:this.samples(instance,r.id),color:['#c93e42','#368364','#2679b6','#da9429','#8467a3'][i%5]})),all=series.flatMap(v=>v.values);
      if(!all.length)return;let low=Math.min(...all.map(v=>v.T_C)),high=Math.max(...all.map(v=>v.T_C)),span=high-low;low-=Math.max(1,span*.06);high+=Math.max(1,span*.06);span=high-low;
      const length=instance.domain[this.data.view.axis==='y'?'height':'width'],offset=this.data.view.offset_m;
      for(let i=0;i<=4;i++){const y=240-i*45;svg.append(svgNode('line',{x1:85,x2:985,y1:y,y2:y,stroke:'#e2e9ed'}),svgNode('text',{x:72,y:y+4,'text-anchor':'end','font-size':11,fill:'#728a98'},fmt(this.temperature(low+span*i/4),2)));}
      series.forEach((s,index)=>{let path='';s.values.forEach((v,i)=>path+=(i&&Math.abs(v.lo-s.values[i-1].hi)<1e-12?' L':'M')+(85+(v.lo+v.hi)/2/length*900)+','+(240-(v.T_C-low)/span*180));if(s.values.length===1){const v=s.values[0];svg.append(svgNode('circle',{cx:85+(v.lo+v.hi)/2/length*900,cy:240-(v.T_C-low)/span*180,r:3,fill:s.color}));}svg.append(svgNode('path',{d:path,fill:'none',stroke:s.color,'stroke-width':2.5}),svgNode('text',{x:85+(index%3)*300,y:23+Math.floor(index/3)*15,'font-size':11,fill:s.color},s.name));});
      svg.append(svgNode('text',{x:85,y:275,'font-size':10,fill:'#728a98'},fmt(offset*1000,2)+' mm'),svgNode('text',{x:985,y:275,'text-anchor':'end','font-size':10,fill:'#728a98'},fmt((length+offset)*1000,2)+' mm'),svgNode('text',{x:50,y:50,'font-size':10,fill:'#728a98'},this.unit()));
      const container=$('selected-summary');container.replaceChildren();const table=document.createElement('table'),head=document.createElement('tr');for(const title of ['Região','Mínima ['+this.unit()+']','Máxima ['+this.unit()+']','Saída do fluido ['+this.unit()+']']){const cell=document.createElement('th');cell.textContent=title;head.append(cell);}table.append(head);
      for(const region of instance.result.region_results){const row=document.createElement('tr');for(const value of [region.name,fmt(this.temperature(region.minimum_C),4),fmt(this.temperature(region.maximum_C),4),region.outlet_C==null?'—':fmt(this.temperature(region.outlet_C),4)]){const cell=document.createElement('td');cell.textContent=value;row.append(cell);}table.append(row);}container.append(table);
    }
    exportCSV(){if(!this.result)return;const rows=[['instancia','nome','componente','fator_geracao','regiao','material','i','j','x_m','y_m','temperatura_C','temperatura_K']];for(const instance of this.result.instances)for(const c of instance.result.cells)rows.push([instance.id,instance.name,instance.component,instance.source_factor,c.region,c.material,c.i,c.j,(c.x0+c.x1)/2,(c.y0+c.y1)/2,c.T_C,c.T_C+273.15]);const quote=value=>'"'+String(value).replaceAll('"','""')+'"';this.download('associacao_temperaturas.csv','\ufeff'+rows.map(row=>row.map(quote).join(',')).join('\r\n')+'\r\n','text/csv;charset=utf-8');}
    async exportPNG(){if(!this.result)return;try{const svg=$('assembly-map').cloneNode(true);svg.setAttribute('xmlns','http://www.w3.org/2000/svg');svg.setAttribute('width',1650);svg.setAttribute('height',975);const url=URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(svg)],{type:'image/svg+xml'}));try{const image=new Image();image.src=url;await image.decode();const canvas=document.createElement('canvas');canvas.width=1650;canvas.height=975;const context=canvas.getContext('2d');context.fillStyle='white';context.fillRect(0,0,1650,975);context.drawImage(image,0,0,1650,975);const blob=await new Promise(resolve=>canvas.toBlob(resolve,'image/png'));this.download('associacao.png',blob,'image/png');}finally{URL.revokeObjectURL(url);}}catch(error){this.toast('Não foi possível exportar a imagem: '+error.message);}}
    async exportReport(){if(!this.result)return;try{const [html,css,js]=await Promise.all(['/assembly','/assembly.css','/assembly.js'].map(async path=>{const r=await fetch(path);if(!r.ok)throw new Error('Falha ao carregar o formato do relatório.');return r.text();}));const close='<'+ '/script>',payload=JSON.stringify({assembly:this.data,result:this.result}).replaceAll('<','\\u003c');let output=html.replace('<link rel="stylesheet" href="/assembly.css">','<style>'+css+'</style>');output=output.replace('<script src="/assembly.js">'+close,'<script id="saved-assembly" type="application/json">'+payload+close+'<script>'+js+close);this.download('associacao_relatorio.html',output,'text/html;charset=utf-8');}catch(error){this.toast(error.message);}}
  }
  new AssemblyApp().init();
})();
