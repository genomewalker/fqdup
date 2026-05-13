#!/usr/bin/env python3
"""Bundle plotly-basic.min.js + HTML template into src/damage_html_assets.hpp."""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC  = os.path.join(HERE, '..', 'src')

def make_raw(content: str, varname: str) -> str:
    for delim in ['DHASSET', 'FQDMG_ASSET', 'FQDMG_HTML', 'RAWDELIM_DH']:
        if f'){delim}' not in content:
            return f'static const char {varname}[] = R"{delim}(\n{content}\n){delim}";\n'
    raise ValueError(f"no safe delimiter found for {varname}")

plotly = open(os.path.join(SRC, 'plotly-basic.min.js')).read()

# HTML up to the point where runtime data is inserted
HTML_PRE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>fqdup damage</title>
<style>
:root{--bg:#fafaf8;--elev:#fff;--ink:#1a1a1a;--muted:#6b6960;--faint:#ccc;
  --c-ct:#d97757;--c-ga:#4a8db8;--c-gt:#9b59b6;--c-ok:#27ae60;--c-bad:#e74c3c;--c-warn:#f39c12;
  font-family:system-ui,sans-serif;font-size:14px;}
body{background:var(--bg);color:var(--ink);margin:0;padding:16px;}
h1{font-size:1.1em;font-weight:600;margin:0 0 4px;}
.sub{color:var(--muted);font-size:.85em;margin-bottom:16px;}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-bottom:16px;}
.card{background:var(--elev);border-radius:6px;padding:12px;box-shadow:0 1px 3px #0001;}
.card h2{font-size:.85em;font-weight:600;color:var(--muted);text-transform:uppercase;
  letter-spacing:.05em;margin:0 0 8px;}
table{border-collapse:collapse;width:100%;font-size:.82em;}
th{text-align:left;color:var(--muted);padding:3px 8px 3px 0;border-bottom:1px solid var(--faint);}
td{padding:3px 8px 3px 0;border-bottom:1px solid #f0efeb;}
.badge{display:inline-block;padding:1px 7px;border-radius:10px;font-size:.78em;font-weight:600;}
.det{background:#d4efdf;color:#1a6035;} .notdet{background:#fdecea;color:#a93226;}
.na{background:#f0f0f0;color:#888;}
.bar-wrap{height:18px;background:#f0f0f0;border-radius:4px;overflow:hidden;margin-top:4px;}
.bar-fill{height:100%;border-radius:4px;transition:width .4s;}
.kv{display:flex;flex-wrap:wrap;gap:4px 16px;font-size:.82em;color:var(--muted);}
.kv span b{color:var(--ink);}
</style>
</head>
<body>
<script>
"""
# After PRE, the caller writes: const D = <JSON>;
# Then POST continues.

HTML_POST = r"""
</script>
<script>
(function(){
  document.title = 'fqdup: ' + (D.sample||'damage');
  var header = '<h1>' + (D.sample||'damage') + '</h1>' +
    '<div class="sub">' + D.n_reads.toLocaleString() + ' reads &nbsp;·&nbsp; ' +
    D.library_type + ' &nbsp;·&nbsp; d5=' + D.d5.toFixed(4) + ' &nbsp; d3=' + D.d3.toFixed(4) + '</div>';
  document.body.insertAdjacentHTML('afterbegin', header);

  // ── Damage curves card ───────────────────────────────────────────────────
  var curvesDiv = document.createElement('div');
  curvesDiv.id = 'curves'; curvesDiv.className = 'card';
  curvesDiv.innerHTML = '<h2>Per-position damage rates</h2>';
  document.body.appendChild(curvesDiv);

  // Build smiley-style damage plot using fitted model curves.
  // Raw empirical rates (t_freq_5prime) mix composition bias across read lengths,
  // so the fitted exponential (bg + d * exp(-lam * p)) gives the clean damage signal.
  var pos = D.pos;  // [1, 2, ..., 15]
  var bg5 = D.bg5 || 0, bg3 = D.bg3 || 0;
  var lam5 = D.lam5 || 0.3, lam3 = D.lam3 || 0.3;

  function modelCurve(bg, dmax, lam) {
    return pos.map(function(p){ return bg + dmax * Math.exp(-lam * (p - 1)); });
  }
  function rawClean(arr) {
    // replace -1 sentinels (low coverage) with null
    return arr ? arr.map(function(v){ return v < 0 ? null : v; }) : null;
  }

  var traces = [];
  // Fitted model lines — clean exponential decay (bg + d * exp(-lam * (p-1)))
  if(D.d5 > 0.001) traces.push({x:pos,y:modelCurve(bg5,D.d5,lam5),name:"5' C→T",mode:'lines+markers',
    line:{color:'#d97757',width:2.5},marker:{size:5}});
  if(D.ga5 && D.d3 > 0.001) traces.push({x:pos,y:modelCurve(bg3,D.d3,lam3),name:"3' G→A",mode:'lines+markers',
    line:{color:'#4a8db8',width:2.5,dash:'dot'},marker:{size:5}});
  if(D.ct3 && D.d3 > 0.001) traces.push({x:pos,y:modelCurve(bg3,D.d3,lam3),name:"3' C→T",mode:'lines+markers',
    line:{color:'#e09070',width:2.5,dash:'dot'},marker:{size:5}});
  if(D.gt5) traces.push({x:pos,y:rawClean(D.gt5),name:"5' G→T",mode:'lines+markers',
    line:{color:'#9b59b6',width:2},marker:{size:4}});
  Plotly.newPlot('curves', traces, {
    height:280, margin:{t:10,b:40,l:50,r:10},
    xaxis:{title:'Position from 5′ end',dtick:1},
    yaxis:{title:'Rate',rangemode:'tozero'},
    legend:{orientation:'h',y:-0.22},
    paper_bgcolor:'#fff',plot_bgcolor:'#fff'
  }, {responsive:true, displayModeBar:false});

  // ── Channels table ──────────────────────────────────────────────────────
  var grid = document.createElement('div');
  grid.className = 'grid'; document.body.appendChild(grid);

  var chanCard = document.createElement('div');
  chanCard.className = 'card'; chanCard.innerHTML = '<h2>Damage channels</h2>';
  grid.appendChild(chanCard);

  var tbl = '<table><tr><th>Ch</th><th>Name</th><th>Status</th><th>Signal</th></tr>';
  (D.channels||[]).forEach(function(c){
    var det = c.detected;
    var appl = c.applicable !== undefined ? c.applicable : true;
    var badge, sig = '';
    if(!appl) badge = '<span class="badge na">N/A</span>';
    else if(det) badge = '<span class="badge det">detected</span>';
    else badge = '<span class="badge notdet">not det.</span>';
    if(c.z_score !== undefined && c.z_score !== null)
      sig = 'z=' + c.z_score.toFixed(1);
    else if(c.d5 !== undefined)
      sig = 'd5=' + c.d5.toFixed(4);
    tbl += '<tr><td><b>'+c.id+'</b></td><td>'+c.name.replace(/_/g,' ')+'</td><td>'+badge+'</td><td>'+sig+'</td></tr>';
  });
  tbl += '</table>';
  chanCard.innerHTML += tbl;

  // ── 8-oxoG panel ───────────────────────────────────────────────────────
  var oxCard = document.createElement('div');
  oxCard.className = 'card'; grid.appendChild(oxCard);
  oxCard.innerHTML = '<h2>8-oxoG estimators</h2>';

  var kvhtml = '<div class="kv">' +
    '<span><b>η̄ (RC log-ratio)</b> ' + (D.oxog_eta_bar !== null ? D.oxog_eta_bar.toFixed(5) : 'N/A') + '</span>' +
    '<span><b>g_hat</b> ' + (D.oxog_g_hat !== null ? D.oxog_g_hat.toFixed(5) : 'N/A') + '</span>' +
    '<span><b>trinuc cosine</b> ' + (D.oxog_cosine !== null ? D.oxog_cosine.toFixed(4) : 'N/A') + '</span>' +
    '<span><b>s_gt</b> ' + D.s_gt.toFixed(5) + '</span>' +
    '</div>';
  oxCard.innerHTML += kvhtml;

  if(D.oxog_g_hat !== null){
    var pct = Math.min(100, D.oxog_g_hat * 500);
    var col = D.oxog_g_hat > 0.01 ? '#9b59b6' : '#ccc';
    oxCard.innerHTML += '<div class="bar-wrap" style="margin-top:10px">' +
      '<div class="bar-fill" style="width:'+pct+'%;background:'+col+'"></div></div>' +
      '<div style="font-size:.75em;color:#999;margin-top:2px">g_hat (scaled ×500)</div>';
  }

  // ── Length-stratified smiley curves ────────────────────────────────────
  if(D.length_bins && D.length_bins.length > 0){
    var lbCard = document.createElement('div');
    lbCard.className = 'card';
    lbCard.style.gridColumn = '1 / -1';
    lbCard.innerHTML = '<h2>Length-stratified damage curves</h2><div id="lenbins"></div>';
    document.body.appendChild(lbCard);

    // Colour ramp: short=warm, long=cool
    var ramp5 = ['#e07050','#c0785a','#8b6090'];
    var ramp3 = ['#6aafd4','#4a8db8','#2a6090'];
    var lbTraces = [];
    D.length_bins.forEach(function(b, i){
      var label = b.lo + '–' + b.hi + ' bp';
      var bg5b  = b.bg5  || 0, bg3b  = b.bg3  || 0;
      var lam5b = b.lam5 || 0.3, lam3b = b.lam3 || 0.3;
      var col5  = ramp5[Math.min(i, ramp5.length-1)];
      var col3  = ramp3[Math.min(i, ramp3.length-1)];
      if(b.d5 > 0.001)
        lbTraces.push({x:pos, y:modelCurve(bg5b, b.d5, lam5b),
          name:"5' "+label, mode:'lines+markers',
          line:{color:col5,width:2,dash:i>0?'dot':'solid'}, marker:{size:4}});
      if(b.d3 > 0.001)
        lbTraces.push({x:pos, y:modelCurve(bg3b, b.d3, lam3b),
          name:"3' "+label, mode:'lines+markers',
          line:{color:col3,width:2,dash:i>0?'dot':'solid'}, marker:{size:4}});
    });
    if(lbTraces.length > 0)
      Plotly.newPlot('lenbins', lbTraces, {
        height:260, margin:{t:10,b:50,l:50,r:10},
        xaxis:{title:'Position from 5′ end',dtick:1},
        yaxis:{title:'Rate',rangemode:'tozero'},
        legend:{orientation:'h',y:-0.28},
        paper_bgcolor:'#fff',plot_bgcolor:'#fff'
      }, {responsive:true, displayModeBar:false});
    else
      lbCard.innerHTML += '<p style="color:#999;font-size:.82em">No bins with detectable damage.</p>';
  }

})();
</script>
</body>
</html>
"""

# Embed Plotly between data-script close and viz-script open so the HTML writer
# only needs to emit: HTML_PRE + "const D={...};\n" + HTML_POST
_VIZ_PREFIX = "\n</script>\n<script>\n"
_viz_body   = HTML_POST[len(_VIZ_PREFIX):] if HTML_POST.startswith(_VIZ_PREFIX) else HTML_POST.lstrip()
HTML_POST_FULL = _VIZ_PREFIX + plotly + "\n</script>\n<script>\n" + _viz_body

out_path = os.path.join(SRC, 'damage_html_assets.hpp')
with open(out_path, 'w') as f:
    f.write('// AUTO-GENERATED by scripts/gen_damage_assets.py — do not edit by hand\n')
    f.write('#pragma once\n\n')
    f.write('namespace fqdup_dmghtml {\n\n')
    f.write(make_raw(HTML_PRE, 'HTML_PRE'))
    f.write('\n')
    f.write(make_raw(HTML_POST_FULL, 'HTML_POST'))
    f.write('\n} // namespace fqdup_dmghtml\n')

print(f'Written: {out_path}')
