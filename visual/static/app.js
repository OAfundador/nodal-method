'use strict';
const $ = id => document.getElementById(id);
const clone = value => JSON.parse(JSON.stringify(value));
const fmt = (value, digits = 2) => Number(value).toLocaleString('pt-BR', {maximumFractionDigits: digits, minimumFractionDigits: digits});
const svgNode = (tag, attributes = {}, text = null) => {
  const element = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [key, value] of Object.entries(attributes)) element.setAttribute(key, value);
  if (text !== null) element.textContent = text;
  return element;
};

class NodalEditor {
  constructor() {
    this.model = null; this.examples = {}; this.selected = null; this.activeMaterial = null;
    this.tool = 'select'; this.mesh = null; this.result = null; this.view = 'geometry';
    this.history = []; this.drag = null; this.zoom = 1; this.busy = false; this.editMaterial = null;
    this.bind();
  }

  async init() {
    try {
      const response = await fetch('/api/examples');
      if (!response.ok) throw new Error('Não foi possível carregar os exemplos.');
      this.examples = await response.json(); this.loadModel(this.examples.composite, false, 'composite');
    } catch (error) { this.toast(error.message); }
  }

  toast(message) {
    $('toast').textContent = message; $('toast').hidden = false;
    clearTimeout(this.toastTimer); this.toastTimer = setTimeout(() => $('toast').hidden = true, 6500);
  }

  setBusy(value, message = '') {
    this.busy = value; document.body.classList.toggle('busy', value);
    for (const id of ['solve-model', 'generate-mesh', 'generate-mesh-side', 'book-tests', 'import-model', 'example-select']) $(id).disabled = value;
    if (message) $('state-label').textContent = message;
  }

  async post(route, payload) {
    const response = await fetch(route, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(payload)});
    const data = await response.json();
    if (!response.ok) throw new Error(data.error || 'Falha no cálculo.');
    return data;
  }

  remember(before) { this.history.push(before); if (this.history.length > 35) this.history.shift(); }
  invalidate() {
    this.mesh = null; this.result = null; this.view = 'geometry';
    $('example-select').value = 'custom';
    if (this.model.benchmark) this.model.benchmark.modified = true;
  }
  change(action) {
    if (!this.model || this.busy) return;
    const before = clone(this.model);
    try {
      action(); this.validateGeometry(); this.remember(before); this.invalidate(); this.render();
    } catch (error) { this.model = before; this.render(); this.toast(error.message); }
  }

  validateGeometry() {
    const {width, height, thickness} = this.model.domain;
    if (![width, height, thickness].every(v => Number.isFinite(v) && v > 0)) throw new Error('As dimensões do domínio devem ser positivas.');
    const regions = this.model.regions;
    for (let i = 0; i < regions.length; i++) {
      const a = regions[i];
      if (![a.x0, a.x1, a.y0, a.y1].every(Number.isFinite) || a.x0 < -1e-12 || a.y0 < -1e-12 || a.x1 > width + 1e-12 || a.y1 > height + 1e-12 || a.x1 <= a.x0 || a.y1 <= a.y0)
        throw new Error('O retângulo deve ter dimensões positivas e ficar dentro do domínio.');
      for (let j = 0; j < i; j++) {
        const b = regions[j];
        if (Math.min(a.x1, b.x1) - Math.max(a.x0, b.x0) > 1e-12 && Math.min(a.y1, b.y1) - Math.max(a.y0, b.y0) > 1e-12)
          throw new Error('As regiões podem encostar, mas não se sobrepor. Encaixe suas bordas na grade.');
      }
    }
  }

  loadModel(data, remember = true, preset = 'custom') {
    const before = this.model ? clone(this.model) : null;
    this.model = clone(data);
    try {
      if (!Array.isArray(this.model.materials) || !Array.isArray(this.model.regions) || !this.model.mesh || !this.model.boundaries)
        throw new Error('O arquivo não tem o formato de um modelo do editor.');
      if (this.model.materials.some(m => !m.id || !m.name || !Number.isFinite(m.k) || m.k <= 0) || this.model.regions.some(r => !this.model.materials.some(m => m.id === r.material)))
        throw new Error('Revise os materiais e suas propriedades no arquivo.');
      this.validateGeometry();
    } catch (error) { this.model = before; throw error; }
    if (before && remember) this.remember(before);
    this.selected = null; this.activeMaterial = this.model.materials[0]?.id || null;
    $('example-select').value = preset;
    this.mesh = null; this.result = null; this.zoom = 1; this.view = 'geometry'; this.setTool('select'); this.render();
  }

  bind() {
    document.querySelectorAll('[data-pane]').forEach(button => button.addEventListener('click', () => {
      document.querySelectorAll('[data-pane]').forEach(b => { const active = b === button; b.setAttribute('aria-selected', String(active)); $('pane-' + b.dataset.pane).hidden = !active; });
    }));
    $('example-select').addEventListener('change', event => {
      if (event.target.value === 'custom') return;
      if (event.target.value === 'empty') {
        const data = clone(this.examples.composite); data.name = 'Novo modelo'; data.regions = []; delete data.benchmark;
        for (const side of Object.keys(data.boundaries)) data.boundaries[side] = {kind: 'adiabatic'};
        this.loadModel(data, true, 'empty');
      } else this.loadModel(this.examples[event.target.value], true, event.target.value);
    });
    $('model-name').addEventListener('change', event => { this.model.name = event.target.value || 'Modelo térmico'; });
    $('new-material').addEventListener('click', () => this.materialDialog());
    $('material-form').addEventListener('submit', event => { event.preventDefault(); this.saveMaterial(); });
    document.querySelectorAll('[data-close]').forEach(button => button.addEventListener('click', () => button.closest('dialog').close()));
    $('draw-material').addEventListener('click', () => this.setTool('draw'));
    $('tool-select').addEventListener('click', () => this.setTool('select'));
    $('tool-draw').addEventListener('click', () => this.setTool('draw'));
    $('undo').addEventListener('click', () => this.undo());
    $('view-geometry').addEventListener('click', () => { this.view = 'geometry'; this.render(); });
    $('view-temperatures').addEventListener('click', () => { if (this.result) { this.view = 'temperatures'; this.render(); } });
    $('zoom-in').addEventListener('click', () => { this.zoom = Math.min(1.8, this.zoom * 1.2); this.draw(); });
    $('zoom-out').addEventListener('click', () => { this.zoom = Math.max(.6, this.zoom / 1.2); this.draw(); });
    $('fit-view').addEventListener('click', () => { this.zoom = 1; this.draw(); });
    for (const id of ['generate-mesh', 'generate-mesh-side']) $(id).addEventListener('click', () => this.compute(false));
    $('solve-model').addEventListener('click', () => this.compute(true));
    $('export-csv').addEventListener('click', () => this.exportCSV());
    $('profile-axis').addEventListener('change', () => this.drawProfile());
    $('save-model').addEventListener('click', () => this.download((this.model.name || 'modelo').replace(/[^\p{L}\p{N}_-]/gu, '_') + '.json', JSON.stringify(this.model, null, 2), 'application/json'));
    $('import-model').addEventListener('click', () => $('import-file').click());
    $('import-file').addEventListener('change', async event => {
      const file = event.target.files[0]; if (!file) return;
      try { if (file.size > 1048576) throw new Error('Use um arquivo JSON de até 1 MB.'); this.loadModel(JSON.parse(await file.text())); this.toast('Modelo aberto.'); }
      catch (error) { this.toast(error.message); } finally { event.target.value = ''; }
    });
    for (const id of ['domain-width', 'domain-height', 'domain-thickness']) $(id).addEventListener('change', () => this.change(() => {
      this.model.domain = {width: Number($('domain-width').value) / 1000, height: Number($('domain-height').value) / 1000, thickness: Number($('domain-thickness').value) / 1000};
    }));
    for (const id of ['mesh-nx', 'mesh-ny', 'mesh-align', 'mesh-h']) $(id).addEventListener('change', () => this.change(() => {
      const nx = Number($('mesh-nx').value), ny = Number($('mesh-ny').value);
      if (![nx, ny].every(n => Number.isInteger(n) && n > 0) || nx * ny > 900) throw new Error('Use divisões inteiras positivas, com até 900 células.');
      this.model.mesh = {nx, ny, align: $('mesh-align').checked, h_internal: $('mesh-h').value === '' ? null : Number($('mesh-h').value)};
    }));
    for (const id of ['region-name', 'region-material', 'region-x0', 'region-y0', 'region-width', 'region-height', 'region-source', 'flow-direction', 'flow-mdot', 'flow-tin'])
      $(id).addEventListener('change', () => this.updateRegion());
    $('delete-region').addEventListener('click', () => this.deleteRegion());
    $('book-tests').addEventListener('click', () => this.verification());
    $('canvas').addEventListener('pointerdown', event => this.pointerDown(event));
    $('canvas').addEventListener('pointermove', event => this.pointerMove(event));
    $('canvas').addEventListener('pointerup', event => this.pointerUp(event));
    $('canvas').addEventListener('pointercancel', () => this.cancelDrag());
    $('canvas').addEventListener('dragover', event => event.preventDefault());
    $('canvas').addEventListener('drop', event => {
      event.preventDefault(); const material = event.dataTransfer.getData('text/plain');
      if (!this.model.materials.some(m => m.id === material)) return;
      const p = this.point(event), w = this.model.domain.width / this.model.mesh.nx * 3, h = this.model.domain.height / this.model.mesh.ny * 3;
      this.change(() => this.addRegion(material, Math.min(p.x, this.model.domain.width - w), Math.min(p.y, this.model.domain.height - h), w, h));
    });
    document.addEventListener('keydown', event => {
      if (/INPUT|TEXTAREA|SELECT/.test(event.target.tagName) || document.querySelector('dialog[open]')) return;
      if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 'z') { event.preventDefault(); this.undo(); }
      if (event.key.toLowerCase() === 'v') this.setTool('select');
      if (event.key.toLowerCase() === 'r') this.setTool('draw');
      if (event.key === 'Delete' || event.key === 'Backspace') { event.preventDefault(); this.deleteRegion(); }
      if (event.key === 'Escape') this.cancelDrag();
    });
  }

  undo() {
    if (!this.history.length || this.busy) return;
    this.model = this.history.pop(); this.selected = null; this.activeMaterial = this.model.materials[0]?.id;
    $('example-select').value = 'custom';
    this.mesh = null; this.result = null; this.view = 'geometry'; this.render(); this.toast('Última alteração desfeita.');
  }

  setTool(tool) {
    if (tool === 'draw' && !this.activeMaterial) { this.toast('Cadastre e selecione um material primeiro.'); return; }
    this.tool = tool; $('tool-select').classList.toggle('active', tool === 'select'); $('tool-draw').classList.toggle('active', tool === 'draw');
    $('canvas').classList.toggle('draw', tool === 'draw');
    $('canvas-status').textContent = tool === 'draw' ? 'Clique e arraste para desenhar. As bordas se encaixam na grade.' : 'Selecione uma região para editar.';
  }

  materialDialog(id = null) {
    const material = this.model.materials.find(m => m.id === id) || {};
    this.editMaterial = id; $('material-dialog-title').textContent = id ? 'Editar material' : 'Novo material';
    for (const key of ['name', 'k', 'rho', 'cp', 'mu']) $('mat-' + key).value = material[key] ?? '';
    $('mat-phase').value = material.phase || 'solid'; $('mat-color').value = material.color || '#7190a7';
    $('mat-k-poly').value = (material.k_poly || []).join(', '); $('material-dialog').showModal(); $('mat-name').focus();
  }

  saveMaterial() {
    const data = {id: this.editMaterial || 'm_' + Date.now().toString(36), name: $('mat-name').value.trim(), phase: $('mat-phase').value, color: $('mat-color').value, k: Number($('mat-k').value)};
    if (!data.name || !Number.isFinite(data.k) || data.k <= 0) { this.toast('Informe um nome e uma condutividade positiva.'); return; }
    for (const key of ['rho', 'cp', 'mu']) if ($('mat-' + key).value !== '') {
      data[key] = Number($('mat-' + key).value); if (!Number.isFinite(data[key]) || data[key] <= 0) { this.toast('As propriedades devem ser positivas.'); return; }
    }
    if ($('mat-k-poly').value.trim()) {
      data.k_poly = $('mat-k-poly').value.split(',').map(Number);
      if (!data.k_poly.length || data.k_poly.length > 6 || data.k_poly.some(v => !Number.isFinite(v))) { this.toast('Informe até 6 coeficientes numéricos separados por vírgulas.'); return; }
    }
    this.change(() => { const i = this.model.materials.findIndex(m => m.id === data.id); if (i < 0) this.model.materials.push(data); else this.model.materials[i] = data; this.activeMaterial = data.id; });
    $('material-dialog').close(); this.toast('Material salvo. Escolha Retângulo para desenhar sua região.');
  }

  addRegion(material, x0, y0, width, height) {
    const id = 'r_' + Date.now().toString(36) + Math.random().toString(36).slice(2, 5);
    this.model.regions.push({id, name: 'Região ' + (this.model.regions.length + 1), material, x0, y0, x1: x0 + width, y1: y0 + height, source: 0});
    this.selected = id;
  }

  updateRegion() {
    if (!this.selected) return;
    this.change(() => {
      const r = this.model.regions.find(r => r.id === this.selected); if (!r) return;
      r.name = $('region-name').value || 'Região'; r.material = $('region-material').value;
      r.x0 = Number($('region-x0').value) / 1000; r.y0 = Number($('region-y0').value) / 1000;
      r.x1 = r.x0 + Number($('region-width').value) / 1000; r.y1 = r.y0 + Number($('region-height').value) / 1000;
      r.source = Number($('region-source').value);
      if (!Number.isFinite(r.source)) throw new Error('Informe uma geração de calor numérica.');
      if (this.model.materials.find(m => m.id === r.material)?.phase === 'fluid')
        r.flow = {direction: $('flow-direction').value, mdot: $('flow-mdot').value === '' ? null : Number($('flow-mdot').value), T_in: $('flow-tin').value === '' ? null : Number($('flow-tin').value)};
    });
  }

  deleteRegion() {
    if (!this.selected) return;
    this.change(() => { this.model.regions = this.model.regions.filter(r => r.id !== this.selected); this.selected = null; });
  }

  render() {
    if (!this.model) return;
    $('model-name').value = this.model.name || 'Modelo térmico';
    $('domain-width').value = +(this.model.domain.width * 1000).toFixed(7); $('domain-height').value = +(this.model.domain.height * 1000).toFixed(7); $('domain-thickness').value = +(this.model.domain.thickness * 1000).toFixed(7);
    $('mesh-nx').value = this.model.mesh.nx; $('mesh-ny').value = this.model.mesh.ny; $('mesh-align').checked = this.model.mesh.align !== false; $('mesh-h').value = this.model.mesh.h_internal ?? '';
    const list = $('material-list'); list.replaceChildren();
    for (const material of this.model.materials) {
      const card = document.createElement('div'); card.className = 'material-card' + (material.id === this.activeMaterial ? ' selected' : ''); card.draggable = true; card.tabIndex = 0; card.setAttribute('role', 'button'); card.setAttribute('aria-label', 'Selecionar material ' + material.name);
      const swatch = document.createElement('span'); swatch.className = 'swatch'; swatch.style.backgroundColor = material.color || '#7190a7';
      const copy = document.createElement('div'); copy.className = 'material-copy'; const name = document.createElement('strong'); name.textContent = material.name;
      const value = document.createElement('small'); value.textContent = (material.k_poly ? 'k(T)' : 'k = ' + fmt(material.k, 2)) + ' W/(m·K) · ' + (material.phase === 'fluid' ? 'Fluido' : 'Sólido'); copy.append(name, value);
      const edit = document.createElement('button'); edit.className = 'edit-material'; edit.textContent = '✎'; edit.setAttribute('aria-label', 'Editar material ' + material.name); edit.addEventListener('click', event => { event.stopPropagation(); this.materialDialog(material.id); });
      const select = () => { this.activeMaterial = material.id; this.render(); };
      card.addEventListener('click', select); card.addEventListener('keydown', event => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); select(); } });
      card.addEventListener('dragstart', event => { this.activeMaterial = material.id; event.dataTransfer.setData('text/plain', material.id); }); card.append(swatch, copy, edit); list.append(card);
    }
    this.renderBoundaries(); this.renderInspector();
    const hasResult = !!this.result;
    $('view-temperatures').disabled = !hasResult; $('view-geometry').classList.toggle('active', this.view === 'geometry'); $('view-temperatures').classList.toggle('active', this.view === 'temperatures');
    $('results-empty').hidden = hasResult; $('results-content').hidden = !hasResult; $('export-csv').disabled = !hasResult;
    $('mesh-badge').textContent = (this.mesh?.nx || this.model.mesh.nx) + ' × ' + (this.mesh?.ny || this.model.mesh.ny) + (this.mesh ? ' · ' + this.mesh.cells.length + ' células' : ' divisões');
    $('state-label').textContent = hasResult ? 'Convergência atingida · campo calculado' : this.mesh ? 'Malha conectada · pronta para resolver' : 'Pronto para gerar a malha';
    $('state-indicator').classList.toggle('ready', !!this.mesh);
    $('canvas-heading').textContent = this.view === 'temperatures' ? 'Campo de temperaturas' : 'Geometria do modelo';
    $('benchmark-note').hidden = !this.model.benchmark;
    if (this.model.benchmark) {
      const probe = this.result?.probes?.[0], reference = this.model.benchmark.target_C;
      $('benchmark-description').textContent = this.model.benchmark.description + (probe ? ' Nesta malha: ' + probe.name + ' = ' + fmt(probe.T_C, 4) + ' °C' + (reference != null ? ' (diferença ' + fmt(probe.T_C - reference, 4) + ' °C).' : '.') : '') + (this.model.benchmark.modified ? ' Modelo editado: a referência vale para os dados originais do exemplo.' : '');
    }
    if (hasResult) {
      $('result-max').textContent = fmt(this.result.maximum_C) + ' °C'; $('result-min').textContent = fmt(this.result.minimum_C) + ' °C';
      $('result-residual').textContent = this.result.residual_W.toExponential(2) + ' W'; $('result-size').textContent = this.result.nodes + ' nós · ' + this.result.cells.length + ' células';
      this.drawProfile();
    }
    this.draw();
  }

  renderBoundaries() {
    const container = $('boundary-fields'); container.replaceChildren();
    const labels = {left: 'Esquerda · x = 0', right: 'Direita · x = L', top: 'Superior · y = H', bottom: 'Inferior · y = 0'};
    for (const [side, label] of Object.entries(labels)) {
      const boundary = this.model.boundaries[side] || {kind: 'adiabatic'};
      const block = document.createElement('div'); block.className = 'boundary-block'; const title = document.createElement('h3'); title.textContent = label; block.append(title);
      const select = document.createElement('select'); select.setAttribute('aria-label', 'Contorno ' + side);
      for (const [value, text] of [['adiabatic', 'Isolado / simetria'], ['temperature', 'Temperatura prescrita'], ['convection', 'Convecção com fluido'], ['radiation', 'Radiação para ambiente'], ['convection_radiation', 'Convecção + radiação'], ['heat_flux', 'Fluxo de calor imposto']]) { const option = document.createElement('option'); option.value = value; option.textContent = text; select.append(option); }
      select.value = boundary.kind; block.append(select);
      select.addEventListener('change', () => this.change(() => { this.model.boundaries[side] = {kind: select.value}; }));
      const fields = boundary.kind === 'temperature' ? [['T', 'Temperatura', '°C']] : boundary.kind === 'heat_flux' ? [['q', 'Fluxo entrando', 'W/m²']] : [];
      if (['convection', 'convection_radiation'].includes(boundary.kind)) fields.push(['h', 'Coeficiente h', 'W/(m²·K)'], ['T_inf', 'Temperatura do fluido', '°C']);
      if (['radiation', 'convection_radiation'].includes(boundary.kind)) fields.push(['emissivity', 'Emissividade', '0 a 1'], ['view_factor', 'Fator de vista', '0 a 1'], ['T_env', 'Temperatura radiativa', '°C']);
      for (const [key, name, unit] of fields) {
        const row = document.createElement('label'); row.textContent = name; const units = document.createElement('span'); units.textContent = unit; const input = document.createElement('input'); input.type = 'number'; input.step = 'any'; input.value = boundary[key] ?? ''; input.placeholder = 'Informe o valor'; input.setAttribute('aria-label', label + ' ' + name);
        input.addEventListener('change', () => this.change(() => { const value = input.value === '' ? null : Number(input.value); if (value !== null && !Number.isFinite(value)) throw new Error('Informe um número válido.'); this.model.boundaries[side][key] = value; })); row.append(units, input); block.append(row);
      }
      container.append(block);
    }
  }

  renderInspector() {
    const region = this.model.regions.find(r => r.id === this.selected);
    $('inspector-empty').hidden = !!region; $('region-fields').hidden = !region;
    if (!region) return;
    $('region-name').value = region.name || region.id; const select = $('region-material'); select.replaceChildren();
    for (const material of this.model.materials) { const option = document.createElement('option'); option.value = material.id; option.textContent = material.name; select.append(option); } select.value = region.material;
    $('region-x0').value = +(region.x0 * 1000).toFixed(7); $('region-y0').value = +(region.y0 * 1000).toFixed(7);
    $('region-width').value = +((region.x1 - region.x0) * 1000).toFixed(7); $('region-height').value = +((region.y1 - region.y0) * 1000).toFixed(7); $('region-source').value = region.source || 0;
    $('region-flow').hidden = this.model.materials.find(m => m.id === region.material)?.phase !== 'fluid';
    $('flow-direction').value = region.flow?.direction || 'y-'; $('flow-mdot').value = region.flow?.mdot ?? ''; $('flow-tin').value = region.flow?.T_in ?? '';
  }

  transform() {
    const {width, height} = this.model.domain, scale = Math.min(830 / width, 450 / height) * this.zoom;
    return {scale, left: (1000 - width * scale) / 2, bottom: (620 + height * scale) / 2, width: width * scale, height: height * scale};
  }

  point(event, snap = true) {
    const p = new DOMPoint(event.clientX, event.clientY).matrixTransform($('canvas').getScreenCTM().inverse());
    const t = this.transform(), d = this.model.domain;
    let x = Math.max(0, Math.min(d.width, (p.x - t.left) / t.scale)), y = Math.max(0, Math.min(d.height, (t.bottom - p.y) / t.scale));
    if (snap && $('snap').checked) { x = Math.round(x / d.width * this.model.mesh.nx) * d.width / this.model.mesh.nx; y = Math.round(y / d.height * this.model.mesh.ny) * d.height / this.model.mesh.ny; }
    return {x: +x.toFixed(12), y: +y.toFixed(12)};
  }

  pointerDown(event) {
    if (!this.model || this.busy || event.button !== 0 || this.view === 'temperatures') return;
    const p = this.point(event), regionId = event.target.closest('[data-region]')?.dataset.region;
    const handle = event.target.dataset.handle;
    if (this.tool === 'draw') this.drag = {type: 'draw', start: p, end: p, before: clone(this.model)};
    else if (regionId) {
      this.selected = regionId; const region = this.model.regions.find(r => r.id === regionId);
      this.drag = {type: handle ? 'resize' : 'move', handle, start: p, original: clone(region), before: clone(this.model)}; this.renderInspector();
    } else { this.selected = null; this.renderInspector(); this.draw(); return; }
    $('canvas').setPointerCapture(event.pointerId); this.draw();
  }

  pointerMove(event) {
    if (!this.model) return;
    const p = this.point(event);
    if (!this.drag) {
      let text = 'x = ' + fmt(p.x * 1000) + ' mm · y = ' + fmt(p.y * 1000) + ' mm';
      if (this.view === 'temperatures' && this.result) { const cell = this.result.cells.find(c => p.x >= c.x0 && p.x <= c.x1 && p.y >= c.y0 && p.y <= c.y1); if (cell) text += ' · ' + fmt(cell.T_C) + ' °C'; }
      $('canvas-status').textContent = text; return;
    }
    const drag = this.drag;
    if (drag.type === 'draw') drag.end = p;
    else {
      const region = this.model.regions.find(r => r.id === this.selected), original = drag.original;
      if (drag.type === 'move') {
        const dx = Math.max(-original.x0, Math.min(this.model.domain.width - original.x1, p.x - drag.start.x));
        const dy = Math.max(-original.y0, Math.min(this.model.domain.height - original.y1, p.y - drag.start.y));
        Object.assign(region, {x0: original.x0 + dx, x1: original.x1 + dx, y0: original.y0 + dy, y1: original.y1 + dy});
      } else {
        region.x0 = drag.handle.includes('w') ? p.x : original.x0; region.x1 = drag.handle.includes('e') ? p.x : original.x1;
        region.y0 = drag.handle.includes('s') ? p.y : original.y0; region.y1 = drag.handle.includes('n') ? p.y : original.y1;
      }
      this.renderInspector();
    }
    this.draw();
  }

  pointerUp(event) {
    if (!this.drag) return;
    const drag = this.drag; this.drag = null;
    try {
      if (drag.type === 'draw') {
        const x0 = Math.min(drag.start.x, drag.end.x), y0 = Math.min(drag.start.y, drag.end.y), w = Math.abs(drag.end.x - drag.start.x), h = Math.abs(drag.end.y - drag.start.y);
        if (w < 1e-10 || h < 1e-10) { this.draw(); return; }
        this.addRegion(this.activeMaterial, x0, y0, w, h);
      }
      this.validateGeometry();
      if (JSON.stringify(this.model) !== JSON.stringify(drag.before)) { this.remember(drag.before); this.invalidate(); }
      if (drag.type === 'draw') this.setTool('select');
    } catch (error) { this.model = drag.before; if (!this.model.regions.some(r => r.id === this.selected)) this.selected = null; this.toast(error.message); }
    if ($('canvas').hasPointerCapture(event.pointerId)) $('canvas').releasePointerCapture(event.pointerId);
    this.render();
  }

  cancelDrag() { if (!this.drag) return; this.model = this.drag.before; this.drag = null; this.selected = null; this.render(); }

  temperatureColor(value) {
    const min = this.result.minimum_C, max = this.result.maximum_C;
    const u = max - min < 1e-10 ? .5 : Math.max(0, Math.min(1, (value - min) / (max - min)));
    const stops = [[43,104,161],[155,201,208],[245,230,160],[213,116,66]], position = u * 3, index = Math.min(2, Math.floor(position)), fraction = position - index;
    return 'rgb(' + stops[index].map((v, k) => Math.round(v + (stops[index + 1][k] - v) * fraction)).join(',') + ')';
  }

  draw() {
    if (!this.model) return;
    const canvas = $('canvas'); canvas.replaceChildren(); const t = this.transform(), domain = this.model.domain;
    const x = value => t.left + value * t.scale, y = value => t.bottom - value * t.scale;
    const defs = svgNode('defs'), pattern = svgNode('pattern', {id: 'dots', width: 16, height: 16, patternUnits: 'userSpaceOnUse'});
    pattern.append(svgNode('circle', {cx: 1, cy: 1, r: .7, fill: '#d5dfe5'})); defs.append(pattern);
    const clip = svgNode('clipPath', {id: 'region-clip'});
    for (const r of this.model.regions) clip.append(svgNode('rect', {x: x(r.x0), y: y(r.y1), width: (r.x1 - r.x0) * t.scale, height: (r.y1 - r.y0) * t.scale}));
    defs.append(clip); canvas.append(defs, svgNode('rect', {width: 1000, height: 620, fill: 'url(#dots)'}));
    canvas.append(svgNode('rect', {x: t.left, y: t.bottom - t.height, width: t.width, height: t.height, fill: '#fff', stroke: '#b6c5cf', 'stroke-width': 1}));
    const heat = this.view === 'temperatures' && this.result;
    for (const region of this.model.regions) {
      const material = this.model.materials.find(m => m.id === region.material);
      canvas.append(svgNode('rect', {x: x(region.x0), y: y(region.y1), width: Math.max(0, (region.x1 - region.x0) * t.scale), height: Math.max(0, (region.y1 - region.y0) * t.scale), fill: material?.color || '#7190a7', 'fill-opacity': .33, stroke: material?.color || '#7190a7', 'stroke-width': 1.5, 'data-region': region.id, class: 'region'}));
    }
    if (heat) for (const cell of this.result.cells) canvas.append(svgNode('rect', {x: x(cell.x0), y: y(cell.y1), width: (cell.x1 - cell.x0) * t.scale + .2, height: (cell.y1 - cell.y0) * t.scale + .2, fill: this.temperatureColor(cell.T_C), 'pointer-events': 'none'}));
    if (this.mesh) {
      const grid = svgNode('g', {'clip-path': 'url(#region-clip)', 'pointer-events': 'none'});
      for (let i = 0; i <= this.mesh.nx; i++) grid.append(svgNode('line', {x1: x(domain.width * i / this.mesh.nx), x2: x(domain.width * i / this.mesh.nx), y1: y(0), y2: y(domain.height), stroke: heat ? '#ffffff45' : '#557b9345', 'stroke-width': .9}));
      for (let j = 0; j <= this.mesh.ny; j++) grid.append(svgNode('line', {x1: x(0), x2: x(domain.width), y1: y(domain.height * j / this.mesh.ny), y2: y(domain.height * j / this.mesh.ny), stroke: heat ? '#ffffff45' : '#557b9345', 'stroke-width': .9}));
      canvas.append(grid);
      if (!heat) for (const contact of this.mesh.contacts) canvas.append(svgNode('line', {x1: x(contact.x0), y1: y(contact.y0), x2: x(contact.x1), y2: y(contact.y1), stroke: '#2e8e79', 'stroke-width': 3, 'pointer-events': 'none'}));
    }
    for (const region of this.model.regions) {
      if (region.x1 <= region.x0 || region.y1 <= region.y0) continue;
      const cx = x((region.x0 + region.x1) / 2), cy = y((region.y0 + region.y1) / 2), w = (region.x1 - region.x0) * t.scale, h = (region.y1 - region.y0) * t.scale;
      if (w > 65 && h > 42 && !heat) {
        canvas.append(svgNode('text', {x: cx, y: cy - 2, 'text-anchor': 'middle', 'font-size': 14, fill: '#385567', 'font-weight': 550, 'pointer-events': 'none'}, region.name || region.id));
        const mat = this.model.materials.find(m => m.id === region.material);
        canvas.append(svgNode('text', {x: cx, y: cy + 19, 'text-anchor': 'middle', 'font-size': 11, fill: '#637e8c', 'pointer-events': 'none'}, mat?.name || 'Material'));
        if (region.source) canvas.append(svgNode('text', {x: cx, y: cy + 38, 'text-anchor': 'middle', 'font-size': 10, fill: '#956c47', 'pointer-events': 'none'}, 'q‴ = ' + Number(region.source).toExponential(2) + ' W/m³'));
      }
      if (region.id === this.selected && !heat) {
        canvas.append(svgNode('rect', {x: x(region.x0), y: y(region.y1), width: w, height: h, fill: 'none', stroke: '#2d7f72', 'stroke-width': 2.4, 'pointer-events': 'none'}));
        for (const [handle, hx, hy] of [['nw',region.x0,region.y1],['ne',region.x1,region.y1],['sw',region.x0,region.y0],['se',region.x1,region.y0]])
          canvas.append(svgNode('rect', {x: x(hx) - 5, y: y(hy) - 5, width: 10, height: 10, rx: 1, fill: 'white', stroke: '#2d7f72', 'stroke-width': 1.6, 'data-region': region.id, 'data-handle': handle, class: 'handle' + (handle === 'ne' || handle === 'sw' ? ' alt' : '')}));
      }
    }
    for (let i = 0; i <= 4; i++) {
      const vx = domain.width * i / 4, vy = domain.height * i / 4;
      canvas.append(svgNode('text', {x: x(vx), y: t.bottom + 28, 'text-anchor': 'middle', 'font-size': 10, fill: '#97a7b3', 'pointer-events': 'none'}, fmt(vx * 1000, 1)));
      canvas.append(svgNode('text', {x: t.left - 20, y: y(vy) + 4, 'text-anchor': 'end', 'font-size': 10, fill: '#97a7b3', 'pointer-events': 'none'}, fmt(vy * 1000, 1)));
    }
    canvas.append(svgNode('text', {x: x(domain.width) + 27, y: t.bottom + 28, 'font-size': 10, fill: '#8f9fac'}, 'x [mm]'));
    canvas.append(svgNode('text', {x: t.left - 23, y: y(domain.height) - 22, 'font-size': 10, fill: '#8f9fac'}, 'y [mm]'));
    for (const probe of (this.model.probes || [])) {
      const computed = this.result?.probes?.find(p => p.name === probe.name);
      canvas.append(svgNode('circle', {cx:x(probe.x),cy:y(probe.y),r:5,fill:'white',stroke:'#374f65','stroke-width':2,'pointer-events':'none'}));
      canvas.append(svgNode('text',{x:x(probe.x)-10,y:y(probe.y)-12,'text-anchor':'end','font-size':12,fill:'#374f65','font-weight':650,'pointer-events':'none'},probe.name+(computed?' · '+fmt(computed.T_C,2)+' °C':'')));
    }
    for (const [side, boundary] of Object.entries(this.model.boundaries)) {
      if (boundary.kind === 'adiabatic') continue;
      const vertical = side === 'left' || side === 'right', bx = side === 'right' ? domain.width : 0, by = side === 'top' ? domain.height : 0;
      const color = boundary.kind === 'temperature' ? (Number(boundary.T) > 50 ? '#cb865b' : '#6d9fba') : '#7caa96';
      canvas.append(svgNode('line', {x1: x(vertical ? bx : 0), x2: x(vertical ? bx : domain.width), y1: y(vertical ? 0 : by), y2: y(vertical ? domain.height : by), stroke: color, 'stroke-width': 4, 'pointer-events': 'none'}));
      if (boundary.kind === 'temperature' && boundary.T != null) canvas.append(svgNode('text', {x: vertical ? x(bx) + (side === 'left' ? 12 : -12) : x(domain.width / 2), y: vertical ? y(domain.height) + 24 : y(by) + (side === 'top' ? 20 : -12), 'text-anchor': vertical ? (side === 'left' ? 'start' : 'end') : 'middle', fill: color, 'font-size': 12, 'font-weight': 600, 'pointer-events': 'none'}, fmt(boundary.T, 0) + ' °C'));
    }
    if (!this.model.regions.length) {
      canvas.append(svgNode('text', {x: 500, y: 295, 'text-anchor': 'middle', fill: '#9aafb8', 'font-size': 19, 'pointer-events': 'none'}, 'Desenhe sua primeira região'));
      canvas.append(svgNode('text', {x: 500, y: 326, 'text-anchor': 'middle', fill: '#a5b8c0', 'font-size': 12, 'pointer-events': 'none'}, 'Escolha um material e use a ferramenta Retângulo.'));
    }
    if (this.drag?.type === 'draw') { const a = this.drag.start, b = this.drag.end; canvas.append(svgNode('rect', {x: x(Math.min(a.x,b.x)), y: y(Math.max(a.y,b.y)), width: Math.abs(a.x-b.x)*t.scale, height: Math.abs(a.y-b.y)*t.scale, fill: '#74b29a44', stroke: '#368f76', 'stroke-width': 2, 'stroke-dasharray': '6 4', 'pointer-events': 'none'})); }
    $('fit-view').textContent = Math.round(this.zoom * 100) + '%'; $('color-legend').hidden = !heat;
    if (heat) { $('legend-min').textContent = fmt(this.result.minimum_C) + ' °C'; $('legend-max').textContent = fmt(this.result.maximum_C) + ' °C'; }
  }

  async compute(solve) {
    if (this.busy) return; if (!this.model.regions.length) { this.toast('Desenhe ao menos uma região.'); return; }
    this.setBusy(true, solve ? 'Resolvendo os balanços térmicos…' : 'Conectando as regiões…');
    try {
      const data = await this.post(solve ? '/api/solve' : '/api/mesh', this.model);
      this.mesh = data;
      if (solve) { this.result = data; this.view = 'temperatures'; this.selected = null; }
      this.render(); this.toast(solve ? 'Solução convergida. Passe o cursor pelo mapa para consultar as temperaturas.' : data.contacts.length + ' faces de células conectam regiões vizinhas.');
    } catch (error) { this.toast(error.message); this.render(); }
    finally { this.setBusy(false); }
  }

  drawProfile() {
    if (!this.result) return;
    const svg = $('profile-chart'); svg.replaceChildren(); const axis = $('profile-axis').value;
    const cells = this.result.cells.filter(c => axis === 'x' ? c.j === Math.floor(this.result.ny / 2) : c.i === Math.floor(this.result.nx / 2)).sort((a,b) => axis === 'x' ? a.i-b.i : a.j-b.j);
    if (!cells.length) { svg.append(svgNode('text', {x:360,y:70,'text-anchor':'middle',fill:'#8c9ca6','font-size':12}, 'O corte central não atravessa nenhuma região.')); return; }
    const low = this.result.minimum_C, high = this.result.maximum_C, span = high-low || 1, length = axis === 'x' ? this.model.domain.width : this.model.domain.height;
    const px = c => 66 + (axis === 'x' ? (c.x0+c.x1)/2 : (c.y0+c.y1)/2)/length*590;
    const py = c => 100-(c.T_C-low)/span*72;
    for (const [fraction,text] of [[0,fmt(low,1)],[.5,fmt((low+high)/2,1)],[1,fmt(high,1)]]) {
      const y=100-fraction*72; svg.append(svgNode('line',{x1:66,x2:656,y1:y,y2:y,stroke:'#e9eef1','stroke-width':1}),svgNode('text',{x:56,y:y+3,'text-anchor':'end','font-size':9,fill:'#97a7b1'},text));
    }
    let path=''; cells.forEach((c,i) => { const prev=cells[i-1]; const connected=prev && (axis==='x' ? c.i===prev.i+1 : c.j===prev.j+1); path+=(connected?' L':' M')+px(c)+','+py(c); });
    svg.append(svgNode('path',{d:path,stroke:'#3c8f82','stroke-width':2,fill:'none'}));
    for(const c of cells) svg.append(svgNode('circle',{cx:px(c),cy:py(c),r:2,fill:'#3c8f82'}));
    svg.append(svgNode('text',{x:66,y:122,'font-size':9,fill:'#97a7b1'},'0'),svgNode('text',{x:656,y:122,'text-anchor':'end','font-size':9,fill:'#97a7b1'},fmt(length*1000,1)+' mm'),svgNode('text',{x:35,y:15,'font-size':9,fill:'#97a7b1'},'°C'));
  }

  download(name, content, type) {
    const url = URL.createObjectURL(new Blob([content], {type})), anchor = document.createElement('a'); anchor.href=url; anchor.download=name; anchor.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  exportCSV() {
    if (!this.result) return;
    const quote = value => '"' + String(value).replaceAll('"','""') + '"';
    const rows = [['regiao','material','x_m','y_m','temperatura_C','temperatura_K']];
    for (const c of this.result.cells) rows.push([c.region,c.material,(c.x0+c.x1)/2,(c.y0+c.y1)/2,c.T_C,c.T_C+273.15]);
    this.download('temperaturas.csv','\ufeff'+rows.map(r=>r.map(quote).join(',')).join('\r\n'),'text/csv;charset=utf-8');
  }

  table(headers, rows) {
    const table=document.createElement('table'),head=document.createElement('thead'),tr=document.createElement('tr');
    for(const text of headers){const th=document.createElement('th');th.textContent=text;tr.append(th);}head.append(tr);table.append(head);
    const body=document.createElement('tbody');for(const row of rows){const tr=document.createElement('tr');for(const value of row){const td=document.createElement('td');td.textContent=value;tr.append(td);}body.append(tr);}table.append(body);return table;
  }
  async articleReport(button) {
    if(this.busy)return;
    button.disabled=true;button.textContent='Gerando curvas localmente…';this.setBusy(true,'Gerando curvas do benchmark…');
    try {
      const data=await this.post('/api/onur',{}), link=document.createElement('a');
      link.href=data.url;link.target='_blank';link.rel='noopener';link.textContent='Abrir curvas e parâmetros do artigo →';button.replaceWith(link);
    }catch(error){button.disabled=false;button.textContent='Tentar gerar as curvas novamente';this.toast(error.message);}
    finally{this.setBusy(false);this.render();}
  }
  async verification() {
    if(this.busy)return;
    $('information-title').textContent='Benchmarks independentes';
    $('information-content').textContent='Calculando localmente…';$('information-dialog').showModal();this.setBusy(true,'Executando os casos de referência…');
    try {
      const data=await this.post('/api/literature',{}), container=$('information-content');container.replaceChildren();
      const description=document.createElement('p');
        description.textContent='Bastos & Andrade (1996), capítulo V, exercícios 1 e 2. Os três chips são modelos concentrados; a esfera usa Euler implícito e Δt = 0,5 s.';container.append(description);
        container.append(this.table(['Caso','Referência [°C]','Calculado [°C]','Erro [°C]','Verificação'],data.results.map(r=>[r.name,fmt(r.reference_C,4),fmt(r.computed_C,4),fmt(r.error_C,4),r.passed?'OK':'Rever'])));
        const note=document.createElement('p');note.textContent='Tolerâncias: 0,05 °C nos chips (arredondamento do livro) e 0,032 °C na esfera. O refinamento temporal é verificado nos testes automatizados.';container.append(note);
        if(data.online?.length){
          const heading=document.createElement('h3');heading.textContent='Referências públicas: NAFEMS e INL/TMAP8';container.append(heading);
          container.append(this.table(['Caso','Referência [°C]','Calculado [°C]','Diferença [°C]','Critério [°C]'],data.online.map(r=>[r.name,fmt(r.reference_C,4),fmt(r.computed_C,4),fmt(r.error_C,4),fmt(r.tolerance_C,4)])));
          const explanation=document.createElement('p');explanation.textContent='T4 usa temperatura da superfície, interpolada entre as faces da malha, com erro espacial reduzido pelo refinamento. T2 verifica a radiação não linear. No caso com geração TMAP8, a métrica é o erro RMS de todo o perfil em relação à solução analítica, com 40 células. As referências e os critérios estão documentados em BENCHMARKS.md.';container.append(explanation);
        }
          if(data.onur){
            const r=data.onur, heading=document.createElement('h3');heading.textContent='Onur & Sivrioglu (1993) · os três fenômenos juntos';container.append(heading);
            const explanation=document.createElement('p');explanation.textContent='Placa com geração, convecção e radiação: Bi=Nr=u=1, Tf/Ti=Te/Ti=0,25. Valores adimensionais θ=T/Ti. A verificação usa a solução analítica permanente e o balanço de energia transiente.';container.append(explanation);
            const guide=document.createElement('p');guide.textContent='Bi compara a convecção na face com a condução interna: valores maiores intensificam a convecção. u mede a geração interna: u=0 sem geração, u=1 geração de referência, u=5 cinco vezes essa geração, com os demais dados iguais. Nr mede a radiação. São parâmetros adimensionais.';container.append(guide);
            container.append(this.table(['Verificação','Resultado'],[['θ permanente · face isolada',fmt(r.theta_insulated,6)],['θ permanente · face resfriada',fmt(r.theta_surface,6)],['Erro máximo do perfil permanente',r.steady_profile_error.toExponential(2)],['Erro de energia por passo (adimensional)',r.transient_balance_error.toExponential(2)],['Critério',r.passed?'OK':'Rever']]));
            const note=document.createElement('p');note.textContent='As figuras do artigo usam Y=(T−Ti)/(Te−Ti), que difere de θ. Há divergência entre alguns patamares impressos e as equações; a auditoria no relatório explicita a leitura aproximada e o balanço exato. As curvas calculadas não são uma digitalização dos gráficos publicados.';container.append(note);
            const button=document.createElement('button');button.className='secondary';button.textContent='Gerar ou abrir curvas do artigo →';button.addEventListener('click',()=>this.articleReport(button));container.append(button);
          }
    }catch(error){$('information-content').textContent=error.message;this.toast(error.message);}
    finally{this.setBusy(false);this.render();}
  }
}

const editor = new NodalEditor();
editor.init();
