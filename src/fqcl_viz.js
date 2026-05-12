// fqcl-viz v5 — cluster genealogy visualizer
// Fonts: system stack only (no IBM Plex Mono or other web fonts)

const MONO = "ui-monospace, 'Cascadia Mono', 'Segoe UI Mono', 'Liberation Mono', monospace";
const TRANSITIONS = new Set(['CT','TC','AG','GA']);
const DAMAGE_MASK = 5;

const COL_HEX = {
  damage:      '#d97757',
  transition:  '#4a8db8',
  transversion:'#b88a3a',
  ins:         '#6b9a5a',
  del:         '#8a6ba0',
  seqerr:      '#a39e8c'
};
const CLASS_ORDER = ['damage','transition','transversion','ins','del','seqerr'];
const CLASS_LABEL = {
  damage: 'Damage (C→T/G→A)', transition: 'Transition', transversion: 'Transversion',
  ins: 'Insertion', del: 'Deletion', seqerr: 'Seq error'
};

// ── edge classification ────────────────────────────────────────────────────
function classifyEdge(e) {
  if (e.damage_like) return 'damage';
  if (TRANSITIONS.has(e.ref + e.alt)) return 'transition';
  return 'transversion';
}
function termZone(pos, L) {
  if (pos < DAMAGE_MASK) return "5'";
  if (pos >= L - DAMAGE_MASK) return "3'";
  return 'int';
}

// ── tooltip ────────────────────────────────────────────────────────────────
const tip = document.getElementById('tooltip');
function showTip(html, ev) {
  tip.innerHTML = html;
  tip.classList.add('show');
  moveTip(ev);
}
function moveTip(ev) {
  const w = tip.offsetWidth, h = tip.offsetHeight;
  tip.style.left = Math.min(window.innerWidth - w - 8, ev.clientX + 14) + 'px';
  tip.style.top  = Math.min(window.innerHeight - h - 8, ev.clientY + 14) + 'px';
}
function hideTip() { tip.classList.remove('show'); }
document.addEventListener('mousemove', ev => {
  if (tip.classList.contains('show')) moveTip(ev);
});

// ── formatters ────────────────────────────────────────────────────────────
const fmtInt = n => n.toLocaleString();
const fmtPct = (n, d=1) => (n * 100).toFixed(d) + '%';

// ── data access ───────────────────────────────────────────────────────────
function getTopClusters() {
  return window.FQCL_DATA.top_members.data.clusters;
}
function getClusterData(id) {
  const r = window.FQCL_DATA.clusters[String(id)];
  return r ? r.data.cluster : null;
}
function getGlobalMeta() {
  const top = getTopClusters();
  if (!top.length) return {};
  const r = window.FQCL_DATA.clusters[String(top[0].cluster_id)];
  if (!r) return {};
  const fqcl = r.fqcl || {};
  return { n_clusters: fqcl.n_clusters, metadata: fqcl.metadata || {}, path: fqcl.path };
}

// ── topbar ────────────────────────────────────────────────────────────────
function renderTopbar() {
  const top = getTopClusters();
  const gm = getGlobalMeta();
  const meta = gm.metadata || {};
  const chips = [['clusters', fmtInt(gm.n_clusters || top.length)], ['showing', fmtInt(top.length)]];
  if (meta.n_input_reads) chips.push(['input reads', fmtInt(meta.n_input_reads)]);
  if (meta.library_type)  chips.push(['library', meta.library_type]);
  if (meta.tool_version)  chips.push(['fqdup', meta.tool_version]);
  document.getElementById('ds-meta').innerHTML = chips.map(([k, v]) =>
    `<div class="kv"><span class="k">${k}</span><span class="v">${v}</span></div>`).join('');
}
document.getElementById('theme-toggle').addEventListener('click', () => {
  document.body.dataset.theme = document.body.dataset.theme === 'dark' ? 'light' : 'dark';
});

// ── sidebar size histogram ─────────────────────────────────────────────────
function renderSizeHist(clusters) {
  const el = document.getElementById('size-hist');
  const W = el.clientWidth || 216, H = 52;
  const pad = { l: 4, r: 4, t: 4, b: 14 };
  const vals = clusters.map(c => c.n_members);
  const maxV = Math.max(...vals, 2), minV = Math.min(...vals, 1);
  const bins = Math.min(24, Math.ceil((W - pad.l - pad.r) / 6));
  const logMin = Math.log10(Math.max(1, minV)), logMax = Math.log10(Math.max(2, maxV));
  const edges = Array.from({length: bins + 1}, (_, i) =>
    Math.pow(10, logMin + (logMax - logMin) * i / bins));
  const hist = new Array(bins).fill(0);
  for (const v of vals) {
    let b = edges.findIndex((e, i) => i < bins && v >= edges[i] && v < edges[i + 1]);
    if (b < 0) b = bins - 1;
    hist[Math.max(0, b)]++;
  }
  const yMax = Math.max(...hist, 1);
  const aw = W - pad.l - pad.r, ah = H - pad.t - pad.b;
  const bw = aw / bins;
  const svg = d3.create('svg').attr('width', W).attr('height', H);
  svg.selectAll('rect').data(hist).join('rect')
    .attr('x', (_, i) => pad.l + i * bw)
    .attr('y', d => pad.t + ah - d * ah / yMax)
    .attr('width', Math.max(1, bw - 0.5))
    .attr('height', d => d * ah / yMax)
    .attr('fill', 'var(--accent)').attr('opacity', 0.65);
  svg.append('text').attr('x', pad.l).attr('y', H - 2)
    .attr('font-size', 9).attr('fill', 'var(--ink-faint)').attr('font-family', MONO)
    .text(fmtInt(minV));
  svg.append('text').attr('x', W - pad.r).attr('y', H - 2)
    .attr('text-anchor', 'end').attr('font-size', 9)
    .attr('fill', 'var(--ink-faint)').attr('font-family', MONO)
    .text(fmtInt(maxV) + ' →');
  el.innerHTML = '';
  el.appendChild(svg.node());
}

// ── nav list (virtual scroll) ─────────────────────────────────────────────
let navAll = [], navFiltered = [], activeId = null;
const ITEM_H = 36;

function renderNav(clusters) {
  navAll = clusters;
  document.getElementById('nav-count').textContent = clusters.length;
  applyNavFilter('');
}
function applyNavFilter(q) {
  navFiltered = q ? navAll.filter(x => String(x.cluster_id).includes(q)) : navAll;
  refreshNavList();
}
function refreshNavList() {
  const list = document.getElementById('nav-list');
  const scrollTop = list.scrollTop, viewH = list.clientHeight || 500;
  const total = navFiltered.length, buf = 8;
  const firstVis = Math.floor(scrollTop / ITEM_H);
  const lastVis  = Math.min(total - 1, Math.ceil((scrollTop + viewH) / ITEM_H));
  const first = Math.max(0, firstVis - buf), last = Math.min(total - 1, lastVis + buf);
  let html = `<div style="height:${first * ITEM_H}px"></div>`;
  for (let i = first; i <= last; i++) {
    const it = navFiltered[i];
    const act = it.cluster_id === activeId ? ' active' : '';
    html += `<div class="nav-item${act}" data-id="${it.cluster_id}">` +
            `<span class="id">${it.cluster_id}</span>` +
            `<span class="n">${fmtInt(it.n_members)}</span></div>`;
  }
  html += `<div style="height:${Math.max(0, total - last - 1) * ITEM_H}px"></div>`;
  list.innerHTML = html;
  list.querySelectorAll('.nav-item').forEach(el =>
    el.addEventListener('click', () => loadCluster(+el.dataset.id)));
}
document.getElementById('nav-filter').addEventListener('input', ev => applyNavFilter(ev.target.value.trim()));
document.getElementById('nav-list').addEventListener('scroll', refreshNavList);

// ── cluster loading ────────────────────────────────────────────────────────
let currentCluster = null;

function loadCluster(id) {
  activeId = id;
  refreshNavList();
  const cl = getClusterData(id);
  if (!cl) {
    document.getElementById('sheet').innerHTML =
      `<div style="padding:24px;color:var(--ink-muted)">Cluster ${id} not embedded in this file.</div>`;
    return;
  }
  currentCluster = cl;
  selectedNodeId = null;
  document.getElementById('cluster-id').textContent = `#${cl.cluster_id}`;
  const scored = cl.edges.filter(e => e.score_evaluated).length;
  const dmg = cl.edges.filter(e => classifyEdge(e) === 'damage').length;
  document.getElementById('cluster-stats').innerHTML =
    `<div class="kv"><span class="k">reads</span><span class="v">${fmtInt(cl.n_members)}</span></div>` +
    `<div class="kv"><span class="k">edges</span><span class="v">${cl.n_edges}</span></div>` +
    `<div class="kv"><span class="k">length</span><span class="v">${cl.parent_len} bp</span></div>` +
    `<div class="kv"><span class="k">scored</span><span class="v">${scored}</span></div>` +
    `<div class="kv"><span class="k">damage</span><span class="v">${dmg}</span></div>`;
  document.getElementById('qc-chips').innerHTML = '';
  renderParentTracks(cl);
  renderHeatmap(cl);
  renderTree(cl);
  renderClassBar(cl);
  renderEvidence(cl);
  renderLedger(cl);
}

// ── parent sequence track ─────────────────────────────────────────────────
function renderParentTracks(cl) {
  const el = document.getElementById('parent-tracks');
  if (!cl.parent_seq) { el.innerHTML = ''; return; }
  const seq = cl.parent_seq, L = seq.length;
  const posCls = new Map(), posEdges = new Map();
  for (const e of cl.edges) {
    const ec = classifyEdge(e);
    if (!posEdges.has(e.pos)) posEdges.set(e.pos, []);
    posEdges.get(e.pos).push(e);
    if (!posCls.has(e.pos) || ec === 'damage') posCls.set(e.pos, ec);
  }
  let ruler = '';
  for (let i = 0; i < L; i++) {
    if (i % 10 === 0) ruler += `<span class="rtick">${i}</span>`;
  }
  let bases = '';
  for (let i = 0; i < L; i++) {
    const b = seq[i], ec = posCls.get(i);
    const inMask = i < DAMAGE_MASK || i >= L - DAMAGE_MASK;
    if (ec) {
      const col = COL_HEX[ec];
      const edges = posEdges.get(i);
      const tt = edges.map(e =>
        `${e.ref}→${e.alt} (${classifyEdge(e)}, ${termZone(e.pos, L)}, n=${e.n_reads})`
      ).join('<br>');
      bases += `<span class="sbase variant" style="--vc:${col}" ` +
               `onmouseenter="showTip('pos ${i}: ${tt}',event)" onmouseleave="hideTip()">${b}</span>`;
    } else if (inMask) {
      bases += `<span class="sbase mask">${b}</span>`;
    } else {
      bases += `<span class="sbase">${b}</span>`;
    }
  }
  el.innerHTML =
    `<div class="seq-wrap">` +
    `<div class="seq-ruler mono">${ruler}</div>` +
    `<div class="seq-bases mono">${bases}</div>` +
    `</div>`;
}

// ── damage heatmap ────────────────────────────────────────────────────────
function renderHeatmap(cl) {
  const el = document.getElementById('heatmap-canvas');
  const L = cl.parent_len;
  if (!L || !cl.edges.length) {
    el.innerHTML = '<div style="padding:8px;color:var(--ink-faint);font-size:11px">no edges</div>';
    return;
  }
  const W = Math.min(el.clientWidth || 700, 1000);
  const cellW = Math.max(2, Math.floor(W / L));
  const actualW = cellW * L, rowH = 11, padT = 2, padB = 14;
  const nodeDepth = new Map([[0, 0]]);
  for (const n of cl.nodes) {
    if (!n.is_parent && n.depth !== undefined) nodeDepth.set(n.id, n.depth);
  }
  const maxDepth = Math.max(...nodeDepth.values(), 0);
  const H = (maxDepth + 1) * rowH + padT + padB;
  const svg = d3.create('svg').attr('width', actualW).attr('height', H).style('display', 'block');
  for (const e of cl.edges) {
    const ec = classifyEdge(e), depth = nodeDepth.get(e.from) ?? 0;
    svg.append('rect')
      .attr('x', e.pos * cellW).attr('y', padT + depth * rowH)
      .attr('width', cellW - 0.5).attr('height', rowH - 1)
      .attr('fill', COL_HEX[ec]).attr('rx', 1).attr('opacity', 0.85)
      .on('mouseenter', ev => showTip(
        `pos ${e.pos} · ${e.ref}→${e.alt} · ${ec} · ${termZone(e.pos, L)} · reads=${e.n_reads}` +
        (e.score_evaluated ? ` · S=${(+e.score).toFixed(2)}` : ''), ev))
      .on('mouseleave', hideTip);
  }
  for (let d = 0; d <= Math.min(maxDepth, 8); d++) {
    svg.append('text').attr('x', 2).attr('y', padT + d * rowH + rowH * 0.75)
      .attr('font-size', 7).attr('font-family', MONO).attr('fill', 'var(--ink-faint)').text(`d${d}`);
  }
  for (let i = 0; i < L; i += 10) {
    svg.append('text').attr('x', i * cellW + 1).attr('y', H - 2)
      .attr('font-size', 7).attr('font-family', MONO).attr('fill', 'var(--ink-faint)').text(i);
  }
  el.innerHTML = '';
  const wrap = document.createElement('div');
  wrap.style.overflowX = 'auto';
  wrap.appendChild(svg.node());
  el.appendChild(wrap);
}

// ── genealogy tree ────────────────────────────────────────────────────────
let selectedNodeId = null, treeCluster = null;

function renderTree(cl) {
  treeCluster = cl;
  const el = document.getElementById('tree-canvas');
  el.innerHTML = '';
  if (!cl.nodes || cl.nodes.length <= 1) {
    el.innerHTML = '<div style="padding:8px;color:var(--ink-faint);font-size:11px">no children</div>';
    document.getElementById('gn-rowcount').textContent = '1 node';
    return;
  }
  const edgeByTo = new Map();
  for (const e of cl.edges) { if (!edgeByTo.has(e.to)) edgeByTo.set(e.to, e); }
  const nodeObjs = new Map();
  for (const n of cl.nodes) nodeObjs.set(n.id, {id: n.id, is_parent: !!n.is_parent, children: []});
  for (const n of cl.nodes) {
    if (!n.is_parent && n.parent !== undefined) {
      const par = nodeObjs.get(n.parent);
      if (par) par.children.push(nodeObjs.get(n.id));
    }
  }
  const hierarchyRoot = d3.hierarchy(nodeObjs.get(0));
  const nNodes = hierarchyRoot.descendants().length;
  document.getElementById('gn-rowcount').textContent = `${nNodes} nodes`;
  const nLeaves = hierarchyRoot.leaves().length;
  const treeDepth = hierarchyRoot.height;
  const mL = 24, mR = 90, mT = 14, mB = 10;
  const spacingY = Math.max(18, Math.min(38, 380 / Math.max(nLeaves, 1)));
  const treeH = Math.max(60, nLeaves * spacingY);
  const treeW = Math.max(280, (treeDepth + 1) * 100);
  const W = treeW + mL + mR, H = treeH + mT + mB;
  const layout = d3.tree().size([treeH, treeW]);
  layout(hierarchyRoot);
  const svg = d3.create('svg').attr('width', W).attr('height', H)
    .style('display', 'block').style('overflow', 'visible');
  const g = svg.append('g').attr('transform', `translate(${mL},${mT})`);

  g.selectAll('path.tree-link').data(hierarchyRoot.links()).join('path')
    .attr('class', 'tree-link')
    .attr('d', d3.linkHorizontal().x(d => d.y).y(d => d.x))
    .attr('fill', 'none')
    .attr('stroke', d => { const e = edgeByTo.get(d.target.data.id); return e ? COL_HEX[classifyEdge(e)] : 'var(--line)'; })
    .attr('stroke-width', d => { const e = edgeByTo.get(d.target.data.id); return e ? Math.max(1, Math.min(4, 0.8 + (e.n_reads || 1) * 0.3)) : 1; })
    .attr('stroke-opacity', 0.6);

  g.selectAll('text.score-badge').data(hierarchyRoot.descendants().slice(1)).join('text')
    .attr('class', 'score-badge')
    .attr('x', d => (d.y + d.parent.y) / 2).attr('y', d => (d.x + d.parent.x) / 2 - 3)
    .attr('text-anchor', 'middle').attr('font-size', 7).attr('font-family', MONO)
    .attr('fill', 'var(--ink-faint)')
    .text(d => { const e = edgeByTo.get(d.data.id); return (e && e.score_evaluated) ? `S=${(+e.score).toFixed(1)}` : ''; });

  const nodeG = g.selectAll('g.tree-node').data(hierarchyRoot.descendants()).join('g')
    .attr('class', 'tree-node')
    .attr('transform', d => `translate(${d.y},${d.x})`)
    .attr('cursor', 'pointer')
    .on('click', (ev, d) => { ev.stopPropagation(); selectNode(d.data.id); })
    .on('mouseenter', (ev, d) => {
      const e = edgeByTo.get(d.data.id);
      if (e) showTip(`node ${d.data.id} · ${e.ref}→${e.alt}@${e.pos} · ${classifyEdge(e)} · ${termZone(e.pos, cl.parent_len)} · n=${e.n_reads}${e.score_evaluated ? ` · S=${(+e.score).toFixed(2)}` : ''}`, ev);
      else showTip('root (representative)', ev);
    })
    .on('mouseleave', hideTip);

  nodeG.append('circle')
    .attr('class', 'tree-node-circle')
    .attr('r', d => d.data.is_parent ? 6 : 4)
    .attr('fill', d => { if (d.data.is_parent) return 'var(--accent)'; const e = edgeByTo.get(d.data.id); return e ? COL_HEX[classifyEdge(e)] : 'var(--ink-faint)'; })
    .attr('stroke', 'var(--bg-elev)').attr('stroke-width', 1.5);

  nodeG.filter(d => !d.children).append('text')
    .attr('x', 9).attr('dominant-baseline', 'middle').attr('font-size', 9).attr('font-family', MONO)
    .attr('fill', 'var(--ink-muted)')
    .text(d => { const e = edgeByTo.get(d.data.id); return e ? `${e.ref}→${e.alt}@${e.pos}` : ''; });

  el.appendChild(svg.node());
}

function selectNode(nodeId) {
  selectedNodeId = nodeId;
  d3.selectAll('#tree-canvas .tree-node-circle')
    .attr('stroke', d => d.data.id === nodeId ? 'var(--ink)' : 'var(--bg-elev)')
    .attr('stroke-width', d => d.data.id === nodeId ? 2.5 : 1.5);
  if (!treeCluster) return;
  const e = treeCluster.edges.find(edge => edge.to === nodeId);
  if (e) highlightLedgerRow(treeCluster.edges.indexOf(e));
}

// ── mutation class bar ────────────────────────────────────────────────────
function renderClassBar(cl) {
  const el = document.getElementById('class-bar');
  const leg = document.getElementById('class-legend');
  const counts = {};
  for (const e of cl.edges) { const c = classifyEdge(e); counts[c] = (counts[c] || 0) + 1; }
  const total = cl.edges.length || 1;
  const W = el.clientWidth || 280, H = 16;
  const svg = d3.create('svg').attr('width', W).attr('height', H);
  let x = 0;
  for (const cls of CLASS_ORDER) {
    const n = counts[cls] || 0;
    if (!n) continue;
    const w = (n / total) * W;
    svg.append('rect').attr('x', x).attr('y', 0).attr('width', Math.max(1, w)).attr('height', H)
      .attr('fill', COL_HEX[cls])
      .on('mouseenter', ev => showTip(`${CLASS_LABEL[cls]}: ${n} (${fmtPct(n/total)})`, ev))
      .on('mouseleave', hideTip);
    x += w;
  }
  el.innerHTML = '';
  el.appendChild(svg.node());
  leg.innerHTML = CLASS_ORDER.filter(c => counts[c]).map(c =>
    `<span class="lg"><span class="sw" style="background:${COL_HEX[c]}"></span>${CLASS_LABEL[c]} (${counts[c]})</span>`
  ).join('');
}

// ── evidence: LR histogram + stats ────────────────────────────────────────
function renderEvidence(cl) {
  const el = document.getElementById('evidence');
  el.innerHTML = '';
  const allEdges = cl.edges, L = cl.parent_len;
  const scored = allEdges.filter(e => e.score_evaluated);

  // column 1: LR score histogram (2fr)
  const histCol = document.createElement('div');
  if (!scored.length) {
    histCol.innerHTML = '<div style="font-size:11px;color:var(--ink-faint)">No scored edges.</div>';
  } else {
    const scores = scored.map(e => +e.score);
    const minS = Math.min(...scores), maxS = Math.max(...scores);
    const mean = scores.reduce((a, b) => a + b, 0) / scores.length;
    const elW = 300, H = 80, pad = {l: 26, r: 8, t: 6, b: 20};
    const domain = [Math.min(minS - 0.5, -3), Math.max(maxS + 0.5, 3)];
    const xScale = d3.scaleLinear().domain(domain).range([pad.l, elW - pad.r]);
    const hist = d3.bin().domain(domain).thresholds(xScale.ticks(20))(scores);
    const yMax = Math.max(...hist.map(b => b.length), 1);
    const yScale = d3.scaleLinear().domain([0, yMax]).range([H - pad.b, pad.t]);
    const svg = d3.create('svg').attr('width', elW).attr('height', H)
      .style('display', 'block').style('overflow', 'visible');
    svg.append('line').attr('x1', xScale(0)).attr('x2', xScale(0))
      .attr('y1', pad.t).attr('y2', H - pad.b)
      .attr('stroke', 'var(--line)').attr('stroke-dasharray', '3,2');
    svg.selectAll('rect.bar').data(hist).join('rect')
      .attr('x', d => xScale(d.x0) + 0.5).attr('y', d => yScale(d.length))
      .attr('width', d => Math.max(0, xScale(d.x1) - xScale(d.x0) - 1))
      .attr('height', d => yScale(0) - yScale(d.length))
      .attr('fill', d => (d.x0 ?? 0) >= 0 ? COL_HEX.damage : COL_HEX.transition)
      .attr('opacity', 0.75)
      .on('mouseenter', (ev, d) => showTip(`S ∈ [${(+d.x0).toFixed(1)}, ${(+d.x1).toFixed(1)}): ${d.length} edges`, ev))
      .on('mouseleave', hideTip);
    svg.append('line').attr('x1', xScale(mean)).attr('x2', xScale(mean))
      .attr('y1', pad.t).attr('y2', H - pad.b)
      .attr('stroke', 'var(--ink-muted)').attr('stroke-dasharray', '2,2');
    svg.append('g').attr('transform', `translate(0,${H - pad.b})`)
      .call(d3.axisBottom(xScale).ticks(5).tickSize(3))
      .call(g => g.select('.domain').remove())
      .call(g => g.selectAll('text').attr('font-size', 8).attr('font-family', MONO).attr('fill', 'var(--ink-muted)'))
      .call(g => g.selectAll('line').attr('stroke', 'var(--line)'));
    svg.append('g').attr('transform', `translate(${pad.l},0)`)
      .call(d3.axisLeft(yScale).ticks(3).tickSize(3))
      .call(g => g.select('.domain').remove())
      .call(g => g.selectAll('text').attr('font-size', 8).attr('font-family', MONO).attr('fill', 'var(--ink-muted)'))
      .call(g => g.selectAll('line').attr('stroke', 'var(--line)'));
    svg.append('text').attr('x', (pad.l + elW - pad.r) / 2).attr('y', H)
      .attr('text-anchor', 'middle').attr('font-size', 8).attr('font-family', MONO)
      .attr('fill', 'var(--ink-faint)').text('LR score');
    histCol.appendChild(svg.node());
    const sub = document.createElement('div');
    sub.style.cssText = `font-family:${MONO};font-size:10px;color:var(--ink-muted);margin-top:4px`;
    sub.textContent = `${scored.length}/${allEdges.length} scored · μ=${mean.toFixed(2)} · [${minS.toFixed(1)}, ${maxS.toFixed(1)}]`;
    histCol.appendChild(sub);
  }
  el.appendChild(histCol);

  // column 2: mutation class breakdown
  const statCol1 = document.createElement('div');
  const classCounts = {};
  for (const e of allEdges) { const c = classifyEdge(e); classCounts[c] = (classCounts[c] || 0) + 1; }
  statCol1.innerHTML = CLASS_ORDER.filter(c => classCounts[c]).map(c =>
    `<div class="ev-stat"><span style="color:${COL_HEX[c]}">${c}</span><span>${classCounts[c]}</span></div>`
  ).join('');
  el.appendChild(statCol1);

  // column 3: terminal zone breakdown
  const statCol2 = document.createElement('div');
  const zones = {"5'": 0, int: 0, "3'": 0};
  for (const e of allEdges) zones[termZone(e.pos, L)]++;
  const dmgEdges = allEdges.filter(e => classifyEdge(e) === 'damage');
  statCol2.innerHTML =
    `<div class="ev-stat"><span>5' edges</span><span>${zones["5'"]}</span></div>` +
    `<div class="ev-stat"><span>3' edges</span><span>${zones["3'"]}</span></div>` +
    `<div class="ev-stat"><span>int edges</span><span>${zones.int}</span></div>` +
    `<div class="ev-stat" style="margin-top:5px"><span>5' damage</span><span>${dmgEdges.filter(e=>e.pos<DAMAGE_MASK).length}</span></div>` +
    `<div class="ev-stat"><span>3' damage</span><span>${dmgEdges.filter(e=>e.pos>=L-DAMAGE_MASK).length}</span></div>`;
  el.appendChild(statCol2);
}

// ── mutation ledger ───────────────────────────────────────────────────────
let ledgerEdges = [];

function renderLedger(cl) {
  const el = document.getElementById('ledger');
  ledgerEdges = cl.edges;
  const L = cl.parent_len;
  const sorted = cl.edges.map((e, i) => ({e, i})).sort((a, b) => a.e.pos - b.e.pos);
  el.innerHTML = sorted.map(({e, i}) => {
    const ec = classifyEdge(e), zone = termZone(e.pos, L);
    const reason = e.damage_like ? 'damage' : (e.score_evaluated ? `S=${(+e.score).toFixed(2)}` : '—');
    return `<div class="ledger-row" data-idx="${i}" onclick="selectLedgerRow(${i})">` +
      `<div><span class="ledger-class-dot" style="background:${COL_HEX[ec]}"></span></div>` +
      `<div>${ec}</div><div class="mono">${e.pos}</div><div class="mono">${e.ref}→${e.alt}</div>` +
      `<div>${zone}</div><div>${reason}</div><div class="mono">${e.n_reads}</div>` +
      `<div class="mono">${fmtPct(e.n_reads / Math.max(cl.n_members, 1))}</div></div>`;
  }).join('');
}

function selectLedgerRow(idx) {
  highlightLedgerRow(idx);
  const e = ledgerEdges[idx];
  if (e) highlightTreeNodeById(e.to);
}
function highlightLedgerRow(idx) {
  document.querySelectorAll('#ledger .ledger-row').forEach(r => r.classList.remove('active'));
  const row = document.querySelector(`#ledger .ledger-row[data-idx="${idx}"]`);
  if (row) { row.classList.add('active'); row.scrollIntoView({behavior: 'smooth', block: 'nearest'}); }
}
function highlightTreeNodeById(nodeId) {
  selectedNodeId = nodeId;
  d3.selectAll('#tree-canvas .tree-node-circle')
    .attr('stroke', d => d.data.id === nodeId ? 'var(--ink)' : 'var(--bg-elev)')
    .attr('stroke-width', d => d.data.id === nodeId ? 2.5 : 1.5);
}

// ── export SVG ────────────────────────────────────────────────────────────
document.getElementById('export-svg').addEventListener('click', () => {
  const svg = document.querySelector('#tree-canvas svg');
  if (!svg) return;
  const blob = new Blob([svg.outerHTML], {type: 'image/svg+xml'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `fqcl-cluster-${activeId || 'tree'}.svg`;
  a.click();
  URL.revokeObjectURL(a.href);
});

// ── color mode controls ───────────────────────────────────────────────────
document.querySelectorAll('.ctrl-group [data-color]').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.ctrl-group [data-color]').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    if (currentCluster) renderTree(currentCluster);
  });
});

// ── panel help buttons ────────────────────────────────────────────────────
const PANEL_HELP = {
  'parent sequence':
    '<b>Parent sequence</b> — the representative read fqdup kept as the cluster output. ' +
    'All other reads in this cluster were identified as PCR duplicates or damage variants of this sequence and discarded. ' +
    'Orange underline = damage position (C→T at 5\' or G→A at 3\'). ' +
    'Blue = transition. Brown = transversion. Shaded zone = terminal damage mask (fqdup ignores C→T here when deciding absorption).',
  'damage fingerprint':
    '<b>Damage fingerprint</b> — each cell is one absorbed read, placed at its mismatch position (x-axis) and ' +
    'tree depth (y-axis, d0=direct children of root, d1=grandchildren…). ' +
    'Color = mutation class. Orange clusters near position 0 or the 3\' end = ancient deamination damage. ' +
    'Blue interior = PCR/sequencing error.',
  'genealogy tree':
    '<b>Genealogy tree</b> — the absorption hierarchy for this cluster. ' +
    'Root (red dot) = the parent sequence. Each child node = a read absorbed as a duplicate. ' +
    'Link color = reason for absorption (orange=damage, blue=transition PCR error, brown=transversion). ' +
    'Link width = number of reads represented. ' +
    'S= label = LR score (log-likelihood ratio; higher = more confident this absorption is a true PCR duplicate). ' +
    'Click a node to highlight its row in the mutation ledger below.',
  'mutation class':
    '<b>Mutation class breakdown</b> — proportion of absorbed edges by type across this cluster. ' +
    'Damage = C→T at 5\' terminal or G→A at 3\' terminal (ancient DNA deamination). ' +
    'Transition = A↔G or C↔T outside the mask zone (likely PCR/seq error). ' +
    'Transversion = any other base change.',
  'evidence':
    '<b>LR score histogram</b> — distribution of log-likelihood ratio scores for all absorbed edges in this cluster. ' +
    'Score > 0 (orange bars, right of dashed line) = absorption more likely a true PCR duplicate than a real variant. ' +
    'Score < 0 (blue bars) = weaker evidence. Dashed vertical line = mean. ' +
    'Right panels: edge counts by mutation class and by terminal zone (5\', interior, 3\').',
  'mutation ledger':
    '<b>Mutation ledger</b> — one row per absorbed edge, sorted by position. ' +
    'Shows: mutation class · position in parent sequence · base change · terminal zone · ' +
    'LR score (S=) · read count · fraction of cluster. ' +
    'Click any row to highlight the corresponding node in the genealogy tree above.'
};

function addPanelHelp() {
  document.querySelectorAll('.panel-head').forEach(head => {
    const title = (head.querySelector('.ph-l1') || {}).textContent || '';
    const key = title.trim().toLowerCase();
    const helpText = PANEL_HELP[key];
    if (!helpText) return;
    const btn = document.createElement('button');
    btn.className = 'help-btn';
    btn.textContent = '?';
    btn.title = 'About this panel';
    btn.addEventListener('click', ev => {
      ev.stopPropagation();
      const existing = document.querySelector('.help-popover');
      if (existing) { existing.remove(); return; }
      const pop = document.createElement('div');
      pop.className = 'help-popover';
      pop.innerHTML = helpText;
      const close = document.createElement('button');
      close.textContent = '✕';
      close.className = 'help-popover-close';
      close.onclick = () => pop.remove();
      pop.appendChild(close);
      document.body.appendChild(pop);
      const r = btn.getBoundingClientRect();
      pop.style.top  = Math.min(window.innerHeight - pop.offsetHeight - 8, r.bottom + 6) + 'px';
      pop.style.left = Math.max(8, Math.min(window.innerWidth - pop.offsetWidth - 8, r.left)) + 'px';
    });
    head.appendChild(btn);
  });
  document.addEventListener('click', () => {
    document.querySelector('.help-popover')?.remove();
  });
}

// ── injected styles for new components ───────────────────────────────────
(function() {
  const s = document.createElement('style');
  s.textContent = `
    .ledger-row.active { background: var(--bg-sunk); outline: 1px solid var(--accent); }
    .seq-wrap { overflow-x: auto; padding: 4px 0; }
    .seq-ruler { font-size: 9px; color: var(--ink-faint); margin-bottom: 2px; white-space: nowrap; }
    .rtick { display: inline-block; width: 80px; text-align: left; }
    .seq-bases { font-size: 12px; line-height: 1.7; white-space: nowrap; letter-spacing: 0.04em; }
    .sbase { display: inline-block; }
    .sbase.mask { background: var(--c-mask); }
    .sbase.variant { background: transparent; border-bottom: 2px solid var(--vc); cursor: help; box-shadow: inset 0 -12px 0 -8px var(--vc); }
    #tree-canvas { overflow-x: auto; }
    .score-badge { pointer-events: none; }
    .help-btn {
      margin-left: auto; width: 16px; height: 16px; border-radius: 50%;
      border: 1px solid var(--line); background: var(--bg-sunk);
      color: var(--ink-muted); font-size: 10px; line-height: 1;
      cursor: pointer; padding: 0; flex-shrink: 0;
    }
    .help-btn:hover { background: var(--accent); color: var(--bg-elev); border-color: var(--accent); }
    .help-popover {
      position: fixed; z-index: 200; max-width: 340px;
      background: var(--bg-elev); border: 1px solid var(--line);
      border-radius: 5px; padding: 12px 14px 10px;
      font-size: 12px; line-height: 1.5; color: var(--ink);
      box-shadow: 0 4px 16px rgba(0,0,0,0.15);
    }
    .help-popover b { color: var(--ink); }
    .help-popover-close {
      position: absolute; top: 6px; right: 8px;
      background: none; border: none; cursor: pointer;
      color: var(--ink-faint); font-size: 12px; padding: 0;
    }
  `;
  document.head.appendChild(s);
})();

// ── init ──────────────────────────────────────────────────────────────────
function init() {
  const clusters = getTopClusters();
  renderTopbar();
  renderSizeHist(clusters);
  renderNav(clusters);
  if (clusters.length > 0) loadCluster(clusters[0].cluster_id);
  addPanelHelp();
}
init();
