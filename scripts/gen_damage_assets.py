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

HTML_PRE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>fqdup damage</title>
<style>
:root{
  --bg:#fafaf8;--elev:#fff;--ink:#1a1a1a;--muted:#6b6960;--faint:#e0dfd8;
  --c-ct:#d97757;--c-ga:#4a8db8;--c-gt:#9b59b6;
  --c-ok:#27ae60;--c-bad:#e74c3c;--c-warn:#f39c12;--c-info:#3498db;
  font-family:system-ui,sans-serif;font-size:14px;
}
*{box-sizing:border-box;}
body{background:var(--bg);color:var(--ink);margin:0;padding:16px 20px;}
h1{font-size:1.15em;font-weight:700;margin:0 0 2px;}
.sub{color:var(--muted);font-size:.85em;margin-bottom:4px;}
.interp{font-size:.85em;color:#444;margin-bottom:14px;font-style:italic;}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-bottom:12px;}
.grid3{display:grid;grid-template-columns:1fr 1fr 1fr;gap:12px;margin-bottom:12px;}
.full{grid-column:1/-1;}
.card{background:var(--elev);border-radius:8px;padding:12px 14px;
      box-shadow:0 1px 4px #0001;border:1px solid var(--faint);}
.card h2{font-size:.78em;font-weight:700;color:var(--muted);text-transform:uppercase;
         letter-spacing:.06em;margin:0 0 10px;}
table{border-collapse:collapse;width:100%;font-size:.82em;}
th{text-align:left;color:var(--muted);padding:3px 8px 3px 0;
   border-bottom:1px solid var(--faint);white-space:nowrap;}
td{padding:3px 8px 3px 0;border-bottom:1px solid #f4f3ef;vertical-align:middle;}
.badge{display:inline-block;padding:1px 7px;border-radius:10px;font-size:.76em;font-weight:600;}
.det  {background:#d4efdf;color:#1a6035;}
.notdet{background:#fdecea;color:#a93226;}
.na   {background:#f0f0f0;color:#888;}
.warn {background:#fef3cd;color:#7d6008;}
.info {background:#d6eaf8;color:#1a5276;}
.pl-excellent{background:#c8f4d8;color:#0d5c2b;}
.pl-good     {background:#d4efdf;color:#1a6035;}
.pl-moderate {background:#fef3cd;color:#7d6008;}
.pl-weak     {background:#fde8c8;color:#8b4513;}
.pl-poor     {background:#fdecea;color:#a93226;}
.pl-none     {background:#f0f0f0;color:#888;}
.kv{display:flex;flex-wrap:wrap;gap:4px 18px;font-size:.82em;color:var(--muted);}
.kv span b{color:var(--ink);}
.bar-wrap{height:10px;background:#f0f0f0;border-radius:4px;overflow:hidden;margin-top:3px;}
.bar-fill{height:100%;border-radius:4px;}
.flags{display:flex;flex-wrap:wrap;gap:4px;margin-bottom:12px;}
.flag-chip{display:inline-block;padding:2px 8px;border-radius:10px;font-size:.76em;
           font-weight:600;background:#fef3cd;color:#7d6008;border:1px solid #f0d080;}
.mono{font-family:monospace;font-size:.88em;}
.dc-yes{color:var(--c-ok);}
.dc-no{color:var(--c-bad);}
</style>
</head>
<body>
<script>
"""

HTML_POST = r"""
</script>
<script>
(function(){
  function jv(v,d){ return (v===null||v===undefined||!isFinite(v))?'N/A':v.toFixed(d||4); }
  function pct(v){ return isFinite(v)?(v*100).toFixed(1)+'%':'N/A'; }
  function badge(cls,text){ return '<span class="badge '+cls+'">'+text+'</span>'; }
  function presClass(label){
    return {excellent:'pl-excellent',good:'pl-good',moderate:'pl-moderate',
            weak:'pl-weak',poor:'pl-poor'}[label]||'pl-none';
  }
  function scoreBar(v,color,max){
    var w=isFinite(v)?Math.min(100,v/(max||1)*100):0;
    return '<div class="bar-wrap"><div class="bar-fill" style="width:'+w+'%;background:'+color+'"></div></div>';
  }
  function appendCard(parent,h2,html,extra){
    var d=document.createElement('div');
    d.className='card'+(extra?' '+extra:'');
    d.innerHTML='<h2>'+h2+'</h2>'+html;
    parent.appendChild(d); return d;
  }

  // ── Header ─────────────────────────────────────────────────────────────────
  document.title='fqdup: '+(D.sample||'damage');
  var pl=(D.pres_label||'').toLowerCase();
  var header='<h1>'+(D.sample||'damage')+'  '+
    (D.pres_label?badge(presClass(pl),D.pres_label):'')+'</h1>'+
    '<div class="sub">'+D.n_reads.toLocaleString()+' reads &nbsp;\xb7&nbsp; '+
    D.library_type+'&nbsp;\xb7&nbsp; d5='+jv(D.d5)+' &nbsp; d3='+jv(D.d3)+'</div>';
  if(D.interpretation)
    header+='<div class="interp">'+D.interpretation+'</div>';
  document.body.insertAdjacentHTML('afterbegin',header);

  if(D.qc_flags&&D.qc_flags.length){
    var fhtml='<div class="flags">';
    D.qc_flags.forEach(function(f){fhtml+='<span class="flag-chip">&#9873; '+f.replace(/_/g,' ')+'</span>';});
    document.body.insertAdjacentHTML('beforeend',fhtml+'</div>');
  }

  // ── Damage curves (full width) ─────────────────────────────────────────────
  var curvesWrap=document.createElement('div');
  curvesWrap.style.marginBottom='12px';
  appendCard(curvesWrap,'Per-position damage rates (fitted model)','<div id="curves"></div>');
  document.body.appendChild(curvesWrap);

  var pos=D.pos, bg5=D.bg5||0, bg3=D.bg3||0, lam5=D.lam5||0.3, lam3=D.lam3||0.3;
  function modelCurve(bg,dmax,lam){
    return pos.map(function(p){return bg+dmax*Math.exp(-lam*(p-1));});
  }
  // 3' curves peak at the rightmost position (3' terminal) — smiley U-shape
  function modelCurveRev(bg,dmax,lam){
    var N=pos[pos.length-1];
    return pos.map(function(p){return bg+dmax*Math.exp(-lam*(N-p));});
  }
  function rawClean(arr){
    return arr?arr.map(function(v){return v<0?null:v;}):null;
  }

  var traces=[];
  if(D.d5>0.001) traces.push({x:pos,y:modelCurve(bg5,D.d5,lam5),name:"5' C→T",
    mode:'lines+markers',line:{color:'#d97757',width:2.5},marker:{size:5}});
  if(D.ga5&&D.d3>0.001) traces.push({x:pos,y:modelCurveRev(bg3,D.d3,lam3),name:"3' G→A",
    mode:'lines+markers',line:{color:'#4a8db8',width:2.5,dash:'dot'},marker:{size:5}});
  if(D.ct3&&D.d3>0.001) traces.push({x:pos,y:modelCurveRev(bg3,D.d3,lam3),name:"3' C→T",
    mode:'lines+markers',line:{color:'#e09070',width:2.5,dash:'dot'},marker:{size:5}});
  if(D.gt5) traces.push({x:pos,y:rawClean(D.gt5),name:"5' G→T",
    mode:'lines+markers',line:{color:'#9b59b6',width:2},marker:{size:4}});
  if(!traces.length) traces.push({x:pos,y:pos.map(function(){return 0;}),
    name:'(no damage)',mode:'lines',line:{color:'#ccc',width:1,dash:'dot'}});
  Plotly.newPlot('curves',traces,{
    height:240,margin:{t:8,b:40,l:50,r:10},
    xaxis:{title:'Position from 5′ end',dtick:1},
    yaxis:{title:'Rate',rangemode:'tozero'},
    legend:{orientation:'h',y:-0.24},
    paper_bgcolor:'#fff',plot_bgcolor:'#fff'
  },{responsive:true,displayModeBar:false});

  // ── Row 1: Channels · 8-oxoG · Library typing ─────────────────────────────
  var row1=document.createElement('div'); row1.className='grid3';
  document.body.appendChild(row1);

  var tbl='<table><tr><th>Ch</th><th>Name</th><th>Status</th><th>Signal</th></tr>';
  (D.channels||[]).forEach(function(c){
    var appl=c.applicable!==undefined?c.applicable:true;
    var b,sig='';
    if(!appl) b=badge('na','N/A');
    else if(c.detected) b=badge('det','detected');
    else b=badge('notdet','not det.');
    if(c.z_score!=null&&isFinite(c.z_score)) sig='z='+jv(c.z_score,1);
    else if(c.d5!=null) sig='d5='+jv(c.d5);
    tbl+='<tr><td><b>'+c.id+'</b></td><td>'+c.name.replace(/_/g,' ')+'</td><td>'+b+'</td><td>'+sig+'</td></tr>';
  });
  appendCard(row1,'Damage channels',tbl+'</table>');

  var oxHtml='<div class="kv">'+
    '<span><b>η̅ (RC log-ratio)</b> '+jv(D.oxog_eta_bar,5)+'</span>'+
    '<span><b>g_hat</b> '+jv(D.oxog_g_hat,5)+'</span>'+
    '<span><b>trinuc cosine</b> '+jv(D.oxog_cosine,4)+'</span>'+
    '<span><b>s_gt</b> '+jv(D.s_gt,5)+'</span>'+
    '</div>';
  if(D.oxog_g_hat!=null&&isFinite(D.oxog_g_hat)){
    var col2=D.oxog_g_hat>0.01?'#9b59b6':'#ccc';
    oxHtml+=scoreBar(D.oxog_g_hat,col2,0.002)+
      '<div style="font-size:.72em;color:#999">g_hat (scaled \xd7500)</div>';
  }
  appendCard(row1,'8-oxoG estimators',oxHtml);

  var libHtml='<div class="kv" style="margin-bottom:8px">'+
    '<span><b>Model</b> '+(D.lib_bic_model||'—')+'</span>'+
    '<span><b>BIC margin</b> '+jv(D.lib_bic_margin,1)+'</span>'+
    '<span><b>p(DS)</b> '+jv(D.lib_p_ds,3)+'</span>'+
    '<span><b>p(SS)</b> '+jv(D.lib_p_ss,3)+'</span>'+
    '</div>';
  if(D.lib_artifact) libHtml+=badge('warn','artifact contamination');
  appendCard(row1,'Library typing',libHtml);

  // ── Row 2: Damage context · Preservation · Depurination ───────────────────
  var row2=document.createElement('div'); row2.className='grid3';
  document.body.appendChild(row2);

  function ctxBar(label,val,color){
    var w=isFinite(val)?Math.min(100,val*100):0;
    return '<div style="display:flex;align-items:center;gap:6px;margin-bottom:5px;font-size:.82em">'+
      '<span style="min-width:110px;color:var(--muted)">'+label+'</span>'+
      '<div class="bar-wrap" style="flex:1"><div class="bar-fill" style="width:'+w+'%;background:'+color+'"></div></div>'+
      '<span style="min-width:38px;text-align:right">'+jv(val,3)+'</span></div>';
  }
  var ctxHtml='';
  if(D.dom_process) ctxHtml+='<div style="margin-bottom:8px">'+badge('info',D.dom_process.replace(/_/g,' '))+'</div>';
  ctxHtml+=ctxBar('deamination',  D.dam_score, '#d97757');
  ctxHtml+=ctxBar('oxidation',    D.ox_score,  '#9b59b6');
  ctxHtml+=ctxBar('fragmentation',D.frag_score,'#4a8db8');
  ctxHtml+=ctxBar('lib. artifact',D.art_score, '#e74c3c');
  ctxHtml+='<div class="kv" style="margin-top:6px">'+
    '<span><b>CpG contrast</b> '+jv(D.cpg_contrast,4)+'</span>'+
    '<span><b>CpG z</b> '+jv(D.cpg_z,1)+'</span></div>';
  appendCard(row2,'Damage context',ctxHtml);

  var presHtml='';
  if(D.pres_score!==undefined&&isFinite(D.pres_score)){
    var sc=D.pres_score;
    var scCol=sc>0.6?'#27ae60':sc>0.3?'#f39c12':'#e74c3c';
    presHtml='<div style="font-size:1.6em;font-weight:700;color:'+scCol+'">'+jv(sc,3)+
      '<span style="font-size:.5em;font-weight:400;color:var(--muted)"> score</span></div>';
    presHtml+=scoreBar(sc,scCol,1)+'<div style="margin-bottom:8px"></div>';
    presHtml+='<div class="kv">'+
      '<span><b>authenticity</b> '+jv(D.auth_eff,3)+'</span>'+
      '<span><b>oxidation</b> '+jv(D.ox_eff,3)+'</span>'+
      '<span><b>QC risk</b> '+jv(D.qcr_eff,3)+'</span>'+
      '</div>';
  }
  appendCard(row2,'Preservation',presHtml);

  var depHtml='<div style="margin-bottom:8px">'+
    (D.depur_detected?badge('warn','detected'):badge('det','not detected'))+'</div>'+
    '<div class="kv">'+
    '<span><b>enrichment 5′</b> '+jv(D.depur_enrich5,5)+'</span>'+
    '<span><b>enrichment 3′</b> '+jv(D.depur_enrich3,5)+'</span>'+
    '<span><b>z-score</b> '+jv(D.depur_z,1)+'</span></div>';
  appendCard(row2,'Depurination',depHtml);

  // ── Library QC (full width) ────────────────────────────────────────────────
  var qcCard=appendCard(document.body,'Library QC &amp; adapter remnants','');
  qcCard.style.marginBottom='12px';
  var qcGrid=document.createElement('div');
  qcGrid.style.cssText='display:grid;grid-template-columns:1fr 1fr;gap:16px';

  var leftDiv=document.createElement('div');
  var leftHtml='';
  if(D.adapter_seq)
    leftHtml+='<div style="margin-bottom:8px">'+badge('warn',D.adapter_seq)+
      ' <span style="font-size:.82em;color:var(--muted)">'+(D.adapter_name||'unknown')+'</span></div>';
  else
    leftHtml+='<div style="margin-bottom:8px">'+badge('det','no adapter stubs detected')+'</div>';
  if(D.stubs5&&D.stubs5.length)
    leftHtml+='<div style="font-size:.82em;margin-bottom:3px"><b>5′ stubs:</b> '+
      D.stubs5.map(function(s){return '<span class="mono">'+s+'</span>';}).join(' ')+'</div>';
  if(D.stubs3&&D.stubs3.length)
    leftHtml+='<div style="font-size:.82em;margin-bottom:3px"><b>3′ stubs:</b> '+
      D.stubs3.map(function(s){return '<span class="mono">'+s+'</span>';}).join(' ')+'</div>';
  leftHtml+='<div class="kv" style="margin-top:8px">'+
    '<span><b>pos0 art. 5′</b> '+(D.pos0_art5?badge('warn','yes'):badge('det','no'))+'</span>'+
    '<span><b>pos0 art. 3′</b> '+(D.pos0_art3?badge('warn','yes'):badge('det','no'))+'</span>'+
    '<span><b>hex entropy 5′</b> '+jv(D.hex_ent5,2)+'</span>'+
    '<span><b>hex entropy int.</b> '+jv(D.hex_ent_int,2)+'</span>'+
    '<span><b>hex JSD</b> '+jv(D.hex_jsd,4)+'</span>'+
    '<span><b>hex shift z</b> '+jv(D.hex_shift_z,1)+'</span>'+
    '<span><b>short-read frac.</b> '+pct(D.short_frac)+'</span>'+
    '</div>';
  leftDiv.innerHTML=leftHtml; qcGrid.appendChild(leftDiv);

  var rightDiv=document.createElement('div');
  if(D.top_hex&&D.top_hex.length){
    var htbl='<div style="font-size:.82em;color:var(--muted);margin-bottom:4px">'+
      'Top enriched 5′ hexamers vs. interior</div>'+
      '<table><tr><th>Hexamer</th><th>log₂ FC</th><th>damage-consistent</th></tr>';
    D.top_hex.forEach(function(h){
      htbl+='<tr><td class="mono">'+h.seq+'</td><td>'+jv(h.log2fc,2)+'</td>'+
        '<td>'+(h.dc?'<span class="dc-yes">yes</span>':'<span class="dc-no">no</span>')+'</td></tr>';
    });
    rightDiv.innerHTML=htbl+'</table>';
  }
  qcGrid.appendChild(rightDiv); qcCard.appendChild(qcGrid);

  // ── Length-stratified smiley curves ────────────────────────────────────────
  if(D.length_bins&&D.length_bins.length){
    var lbCard=appendCard(document.body,'Length-stratified damage curves','<div id="lenbins"></div>');
    lbCard.style.marginBottom='12px';
    var ramp5=['#e07050','#c0785a','#8b6090'];
    var ramp3=['#6aafd4','#4a8db8','#2a6090'];
    var lbTraces=[];
    D.length_bins.forEach(function(b,i){
      var label=b.lo+'–'+b.hi+' bp';
      var bg5b=b.bg5||0,bg3b=b.bg3||0,lam5b=b.lam5||0.3,lam3b=b.lam3||0.3;
      var c5=ramp5[Math.min(i,ramp5.length-1)],c3=ramp3[Math.min(i,ramp3.length-1)];
      if(b.d5>0.001) lbTraces.push({x:pos,y:modelCurve(bg5b,b.d5,lam5b),
        name:"5' "+label,mode:'lines+markers',
        line:{color:c5,width:2,dash:i?'dot':'solid'},marker:{size:4}});
      if(b.d3>0.001) lbTraces.push({x:pos,y:modelCurveRev(bg3b,b.d3,lam3b),
        name:"3' "+label,mode:'lines+markers',
        line:{color:c3,width:2,dash:i?'dot':'solid'},marker:{size:4}});
    });
    if(lbTraces.length)
      Plotly.newPlot('lenbins',lbTraces,{
        height:240,margin:{t:8,b:50,l:50,r:10},
        xaxis:{title:'Position from 5′ end',dtick:1},
        yaxis:{title:'Rate',rangemode:'tozero'},
        legend:{orientation:'h',y:-0.3},
        paper_bgcolor:'#fff',plot_bgcolor:'#fff'
      },{responsive:true,displayModeBar:false});
    else
      lbCard.insertAdjacentHTML('beforeend',
        '<p style="color:#999;font-size:.82em">No bins with detectable damage.</p>');
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
