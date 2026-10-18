'use strict';

// Operações geométricas compartilhadas pelo desenho e pelos testes sem navegador.
const GeometryTools = (() => {
  const overlap = (a0, a1, b0, b1) => Math.min(a1, b1) - Math.max(a0, b0) > 1e-12;
  function snap(value, targets, tolerance, step = 0) {
    const target = targets.reduce((best, edge) => Math.abs(edge - value) < Math.abs(best - value) ? edge : best, Infinity);
    if (Math.abs(target - value) <= tolerance) return target;
    return step ? Math.round(value / step) * step : value;
  }
  function stretch(region, regions, direction = 'nearest') {
    const candidates = [];
    const blocked = new Set();
    for (const other of regions) {
      if (other.id === region.id) continue;
      if (overlap(region.y0, region.y1, other.y0, other.y1)) {
        if (Math.abs(other.x1 - region.x0) <= 1e-12) blocked.add('left');
        if (Math.abs(other.x0 - region.x1) <= 1e-12) blocked.add('right');
        if (other.x1 < region.x0 - 1e-12) candidates.push({side: 'left', key: 'x0', value: other.x1, gap: region.x0 - other.x1});
        if (other.x0 > region.x1 + 1e-12) candidates.push({side: 'right', key: 'x1', value: other.x0, gap: other.x0 - region.x1});
      }
      if (overlap(region.x0, region.x1, other.x0, other.x1)) {
        if (Math.abs(other.y1 - region.y0) <= 1e-12) blocked.add('bottom');
        if (Math.abs(other.y0 - region.y1) <= 1e-12) blocked.add('top');
        if (other.y1 < region.y0 - 1e-12) candidates.push({side: 'bottom', key: 'y0', value: other.y1, gap: region.y0 - other.y1});
        if (other.y0 > region.y1 + 1e-12) candidates.push({side: 'top', key: 'y1', value: other.y0, gap: other.y0 - region.y1});
      }
    }
    const target = candidates.filter(c => !blocked.has(c.side) && (direction === 'nearest' || c.side === direction)).sort((a, b) => a.gap - b.gap)[0];
    if (!target) throw new Error('Não há uma região separada por uma folga nessa direção. Escolha outra borda; faces já encostadas não precisam de ajuste.');
    return {...region, [target.key]: target.value};
  }
  function move(original, dx, dy, regions, domain, tolerance, grid = null) {
    let x0 = original.x0 + dx, y0 = original.y0 + dy;
    const w = original.x1 - original.x0, h = original.y1 - original.y0;
    const neighbors = regions.filter(r => r.id !== original.id);
    const xTargets = [0, domain.width - w];
    for (const r of neighbors) if (overlap(y0, y0 + h, r.y0, r.y1)) xTargets.push(r.x0 - w, r.x1);
    x0 = snap(x0, xTargets, tolerance.x, grid?.x || 0);
    const yTargets = [0, domain.height - h];
    for (const r of neighbors) if (overlap(x0, x0 + w, r.x0, r.x1)) yTargets.push(r.y0 - h, r.y1);
    y0 = snap(y0, yTargets, tolerance.y, grid?.y || 0);
    x0 = Math.max(0, Math.min(domain.width - w, x0));
    y0 = Math.max(0, Math.min(domain.height - h, y0));
    return {...original, x0, x1: x0 + w, y0, y1: y0 + h};
  }
  return {snap, stretch, move};
})();
if (typeof module !== 'undefined') module.exports = GeometryTools;
