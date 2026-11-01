# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
# Regenerates the three pages from launchers.json:
#   python3 build_viz.py launchers.json .
import json, sys
src, outdir = sys.argv[1], sys.argv[2]
d = json.load(open(src))
data = "const LD=" + json.dumps(d, separators=(',', ':')) + ";"
HEAD = '''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>%s</title>
<!-- SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial -->
<!-- Copyright (C) 2026 Franklin Amador. Data: go run ./examples/ring_towers -phasors (kronspice transient simulation of each launcher alone). -->
<style>
:root{--bg:#faf9f5;--card:#f1efe8;--text:#1f1e1d;--muted:#5f5e5a;--border:rgba(0,0,0,.15);--accent:#185FA5;color-scheme:light}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#1f1e1d;--card:#2c2c2a;--text:#f1efe8;--muted:#b4b2a9;--border:rgba(255,255,255,.2);--accent:#85B7EB;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#1f1e1d;--card:#2c2c2a;--text:#f1efe8;--muted:#b4b2a9;--border:rgba(255,255,255,.2);--accent:#85B7EB;color-scheme:dark}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);font:15px/1.6 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:760px;margin:0 auto;padding:20px 16px 32px}h1{font-size:20px;font-weight:500;margin:0 0 4px}p.sub{color:var(--muted);margin:0 0 16px;font-size:14px}
button{font:inherit;font-size:13px;color:var(--text);background:transparent;border:1px solid var(--border);border-radius:8px;padding:6px 12px;cursor:pointer}
button.on{border-color:var(--accent);color:var(--accent)}input[type=range]{width:100%%}
.row{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:0 0 10px;font-size:13px;color:var(--muted)}
.grid{display:grid;grid-template-columns:150px minmax(0,1fr) 64px;gap:8px 12px;align-items:center;font-size:13px;color:var(--muted);margin:0 0 10px}
.grid output{color:var(--text);font-weight:500}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-top:12px}
.card{background:var(--card);border-radius:8px;padding:10px 14px}.card div:first-child{font-size:13px;color:var(--muted)}.card div:last-child{font-size:17px;font-weight:500}
.sw{display:inline-block;width:14px;height:3px;vertical-align:middle;margin-right:4px}canvas{width:100%%;height:auto;display:block}
.note{font-size:12px;color:var(--muted);margin-top:8px}
</style></head><body><main>
'''
COMMON = data + r'''
const N=96,D2R=Math.PI/180,AMP=0.01,$=id=>document.getElementById(id);
const K={A:"A east",B:"B west",C:"C east",D:"D west"};
function rot(v,d){const c=Math.cos(d),s=Math.sin(d);return [v[0]*c+v[1]*s,v[1]*c-v[0]*s];}
function field(o,names){const out=Array.from({length:N},()=>[0,0]);(names||Object.keys(o)).forEach(n=>LD[K[n]].V.forEach((v,k)=>{const r=rot(v,o[n]*D2R);out[k][0]+=r[0];out[k][1]+=r[1];}));return out;}
function power(o,tot){let P=0;Object.keys(o).forEach(n=>LD[K[n]].rad.forEach(r=>{const f=(r.phase+o[n])*D2R;P+=0.5*(tot[r.node][0]*AMP*Math.sin(f)+tot[r.node][1]*AMP*Math.cos(f));}));return P;}
const val=(a,k,t)=>a[k][0]*Math.cos(t)+a[k][1]*Math.sin(t),env=(a,k)=>Math.hypot(a[k][0],a[k][1]);
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const cv=$('cv'),cx=cv.getContext('2d');let ph=0,run=true,last=null;
function frame(VS,curves,envs,ticks){const tp=css('--text'),ts=css('--muted'),bd=css('--border');cx.clearRect(0,0,680,400);
const ox=180,oy=200,R=90,S=55;cx.strokeStyle=bd;cx.lineWidth=1;cx.setLineDash([3,4]);cx.beginPath();cx.arc(ox,oy,R,0,2*Math.PI);cx.stroke();cx.setLineDash([]);
curves.forEach(([f,col,w])=>{cx.strokeStyle=col;cx.lineWidth=w;cx.beginPath();for(let k=0;k<=N;k++){const kk=k%N,a=kk*2*Math.PI/N,r=R+S*f(kk)/VS;const x=ox+r*Math.cos(a),y=oy-r*Math.sin(a);k?cx.lineTo(x,y):cx.moveTo(x,y);}cx.stroke();});
(ticks||[]).forEach(d=>{const a=d*D2R;cx.strokeStyle=tp;cx.lineWidth=2;[-0.03,0.03].forEach(e=>{cx.beginPath();cx.moveTo(ox+(R-10)*Math.cos(a+e),oy-(R-10)*Math.sin(a+e));cx.lineTo(ox+(R+10)*Math.cos(a+e),oy-(R+10)*Math.sin(a+e));cx.stroke();});});
cx.font='13px sans-serif';cx.textAlign='center';cx.textBaseline='middle';
[["A",0,'#7F77DD'],["C",90,'#D85A30'],["B",180,'#7F77DD'],["D",270,'#D85A30']].forEach(([l,d,col])=>{const a=d*D2R;cx.fillStyle=col;cx.beginPath();cx.arc(ox+R*Math.cos(a),oy-R*Math.sin(a),4,0,2*Math.PI);cx.fill();cx.fillStyle=tp;cx.fillText(l+' '+d+'°',ox+(R+S+16)*Math.cos(a),oy-(R+S+16)*Math.sin(a));});
const x0=372,W=294,y0=200,H=150;cx.strokeStyle=bd;cx.lineWidth=0.5;cx.beginPath();cx.moveTo(x0,y0);cx.lineTo(x0+W,y0);cx.stroke();
[0,90,180,270,360].forEach(d=>{const x=x0+W*d/360;cx.beginPath();cx.moveTo(x,y0-H);cx.lineTo(x,y0+H);cx.stroke();cx.fillStyle=ts;cx.font='11px sans-serif';cx.fillText(d+'°',x,y0+H+12);});
cx.fillStyle=ts;cx.textAlign='left';cx.fillText('+'+VS+' V',x0+2,y0-H+8);cx.fillText('−'+VS+' V',x0+2,y0+H-8);cx.textAlign='center';
const lin=(f,col,w,dash)=>{cx.strokeStyle=col;cx.lineWidth=w;cx.setLineDash(dash||[]);cx.beginPath();for(let k=0;k<=N;k++){const kk=k%N,x=x0+W*k/N,y=y0-H*f(kk)/VS;k?cx.lineTo(x,y):cx.moveTo(x,y);}cx.stroke();cx.setLineDash([]);};
(envs||[]).forEach(f=>lin(f,'#888780',1,[4,3]));curves.forEach(([f,col,w])=>lin(f,col,w));}
function controls(){$('play').onclick=()=>{run=!run;$('play').textContent=run?'Pause':'Play';};}
function loop(draw){const step=ts=>{if(last!==null&&run)ph+=(ts-last)/1000*Math.PI*parseFloat($('spd').value);last=ts;draw();requestAnimationFrame(step);};requestAnimationFrame(step);}
'''
PLAY = '<button id="play">Pause</button><label for="spd">Speed</label><input type="range" id="spd" min="0.1" max="2" step="0.1" value="0.6" style="width:110px">'
LEG = lambda items: '<div class="row">' + ''.join(f'<span><span class="sw" style="background:{c}"></span>{t}</span>' for t,c in items) + '<span><span class="sw" style="border-top:1.5px dashed #888780;height:0"></span>envelope</span></div>'

files = {}
# 1. four towers, preset states
files['four_towers.html'] = HEAD % 'Four-tower ring' + '''<h1>Four towers on the ring</h1>
<p class="sub">A (0°), C (90°), B (180°), D (270°); each pair an east and a west launcher. Voltage around the ring from the kronspice simulation, 4 wavelengths around, 0.2 Np loss per lap.</p>
<div class="row" id="btns"></div>
<div class="row">''' + PLAY + '''<label><input type="checkbox" id="comp" checked> show each pair</label></div>
''' + LEG([("total",'#1D9E75'),("A–B pair",'#7F77DD'),("C–D pair",'#D85A30')]) + '''
<canvas id="cv" width="680" height="400"></canvas>
<div class="cards"><div class="card"><div>|V| 0° / 90° / 180° / 270°</div><div id="o-q">–</div></div><div class="card"><div>Towers' power (= ring loss)</div><div id="o-p">–</div></div></div>
<script>''' + COMMON + r'''
const ST={"A–B standing":{A:0,B:0},"C–D standing":{C:0,D:0},"both, synchronized":{A:0,B:0,C:0,D:0},"both, C–D antiphase":{A:0,B:0,C:180,D:180},"both, C–D shifted 1/8 wave":{A:0,B:0,C:-45,D:45}};
let cur="both, synchronized";const bd=$('btns');
Object.keys(ST).forEach(n=>{const b=document.createElement('button');b.textContent=n;b.onclick=()=>{cur=n;sel();};b.dataset.n=n;bd.appendChild(b);});
let F,AB,CD;function sel(){[...bd.children].forEach(b=>b.classList.toggle('on',b.dataset.n===cur));const o=ST[cur];F=field(o);
AB=field({A:o.A||0,B:o.B||0});CD=field({C:o.C||0,D:o.D||0});
$('o-q').textContent=[0,24,48,72].map(k=>env(F,k).toFixed(2)).join(' / ')+' V';$('o-p').textContent=(power(o,F)*1e3).toFixed(1)+' mW';}
sel();controls();
loop(()=>{const both=cur.startsWith('both')&&$('comp').checked;const c=[];if(both){c.push([k=>val(AB,k,ph),'#7F77DD',1.2],[k=>val(CD,k,ph),'#D85A30',1.2]);}c.push([k=>val(F,k,ph),'#1D9E75',2.5]);
frame(7,c,[k=>env(F,k),k=>-env(F,k)]);});
</script></main></body></html>
'''
# 2. phase control
files['phase_control.html'] = HEAD % 'Ring phase control' + '''<h1>Phase control of the two tower pairs</h1>
<p class="sub">Timing and position phase of the C–D pair against A–B. The ring is linear, so the total is the sum of each launcher's simulated response, rotated by your phases.</p>
<div class="row"><button id="p-sync">Synchronized</button><button id="p-anti">C–D antiphase</button><button id="p-shift">C–D shifted 1/8 wave</button><button id="p-quad">C–D quarter cycle late</button></div>
<div class="grid"><label for="tim">C–D timing phase</label><input type="range" id="tim" min="-180" max="180" step="5" value="0"><output id="tim-o"></output>
<label for="pos">C–D position phase</label><input type="range" id="pos" min="-180" max="180" step="5" value="0"><output id="pos-o"></output>
<label for="abt">A–B timing phase</label><input type="range" id="abt" min="-180" max="180" step="5" value="0"><output id="abt-o"></output></div>
<div class="row"><label><input type="checkbox" id="onab" checked> A–B on</label><label><input type="checkbox" id="oncd" checked> C–D on</label><label><input type="checkbox" id="comp" checked> show each pair</label>''' + PLAY + '''</div>
''' + LEG([("total",'#1D9E75'),("A–B pair",'#7F77DD'),("C–D pair",'#D85A30')]) + '''
<canvas id="cv" width="680" height="400"></canvas>
<div class="cards"><div class="card"><div>|V| 0° / 90° / 180° / 270°</div><div id="m-q">–</div></div><div class="card"><div>Peak anywhere</div><div id="m-pk">–</div></div>
<div class="card"><div>Towers' power (= ring loss)</div><div id="m-p">–</div></div><div class="card"><div>Quarter points in step</div><div id="m-s">–</div></div></div>
<div class="note" id="m-l"></div>
<script>''' + COMMON + r'''
const syncP=(()=>{const o={A:0,B:0,C:0,D:0};return power(o,field(o));})();let C={};
function offs(){const t=+$('tim').value,p=+$('pos').value,ab=+$('abt').value;return {A:ab,B:ab,C:t-p,D:t+p};}
function recompute(){const o=offs(),act={};if($('onab').checked){act.A=o.A;act.B=o.B;}if($('oncd').checked){act.C=o.C;act.D=o.D;}
const tot=field(act);C={tot,ab:field({A:o.A,B:o.B}),cd:field({C:o.C,D:o.D}),act};const P=Object.keys(act).length?power(act,tot):0;
['tim','pos','abt'].forEach(i=>$(i+'-o').textContent=$(i).value+'°');
$('m-q').textContent=[0,24,48,72].map(k=>env(tot,k).toFixed(1)).join(' / ')+' V';let pk=0;for(let k=0;k<N;k++)pk=Math.max(pk,env(tot,k));$('m-pk').textContent=pk.toFixed(2)+' V';
$('m-p').textContent=(P*1e3).toFixed(1)+' mW ('+(P/syncP*100).toFixed(0)+'% of sync)';
let bt=0,bk=0;for(let i=0;i<64;i++){const t=i*Math.PI/32,v=Math.abs(val(tot,0,t));if(v>bk){bk=v;bt=t;}}const vs=[0,24,48,72].map(k=>val(tot,k,bt));
const lo=Math.min(...vs.map(Math.abs)),hi=Math.max(...vs.map(Math.abs));$('m-s').textContent=(hi>0.2&&vs.every(v=>v*vs[0]>0))?('yes, spread '+Math.round(100*(hi-lo)/hi)+'%'):'no';
$('m-l').textContent='Launcher phase offsets: A east '+o.A+'°, B west '+o.B+'°, C east '+o.C+'°, D west '+o.D+'°. C–D pattern moved '+(-(+$('pos').value)/4).toFixed(1)+'° around the ring.';}
['tim','pos','abt','onab','oncd'].forEach(id=>$(id).addEventListener('input',recompute));
const set=(t,p)=>{$('tim').value=t;$('pos').value=p;$('abt').value=0;$('onab').checked=true;$('oncd').checked=true;recompute();};
$('p-sync').onclick=()=>set(0,0);$('p-anti').onclick=()=>set(180,0);$('p-shift').onclick=()=>set(0,45);$('p-quad').onclick=()=>set(-90,0);
recompute();controls();
loop(()=>{const c=[],comp=$('comp').checked;if(comp&&'A' in C.act)c.push([k=>val(C.ab,k,ph),'#7F77DD',1.2]);if(comp&&'C' in C.act)c.push([k=>val(C.cd,k,ph),'#D85A30',1.2]);
c.push([k=>val(C.tot,k,ph),'#1D9E75',2.5]);frame(7,c,[k=>env(C.tot,k),k=>-env(C.tot,k)]);});
</script></main></body></html>
'''
# 3. dc bias
files['dc_bias.html'] = HEAD % 'Ring DC bias' + '''<h1>DC bias on each tower pair</h1>
<p class="sub">A continuous ring is one conductor at DC; a sectored ring (10 µF blocking capacitors at 45°, 135°, 225°, 315°) lets each pair hold its own bias. Feeds 1 MΩ. The RF wave rides unchanged on the DC level (checked in kronspice, TestBias).</p>
<div class="row"><button id="m-cont">Continuous ring</button><button id="m-sec">Sectored ring (blocking caps)</button></div>
<div class="grid"><label for="vab">A–B DC bias</label><input type="range" id="vab" min="-20" max="20" step="1" value="10"><output id="vab-o"></output>
<label for="vcd">C–D DC bias</label><input type="range" id="vcd" min="-20" max="20" step="1" value="-5"><output id="vcd-o"></output>
<label for="tim">C–D timing phase</label><input type="range" id="tim" min="-180" max="180" step="5" value="0"><output id="tim-o"></output></div>
<div class="row">''' + PLAY + '''</div>
''' + LEG([("DC + RF",'#1D9E75'),("DC level",'#BA7517')]) + '''
<canvas id="cv" width="680" height="400"></canvas>
<div class="cards"><div class="card"><div>DC at 0° / 90° / 180° / 270°</div><div id="o-dc">–</div></div><div class="card"><div>Peak gap voltage</div><div id="o-pk">–</div></div>
<div class="card"><div>DC current between pairs</div><div id="o-i">–</div></div><div class="card"><div>Steady DC loss in feeds</div><div id="o-p">–</div></div></div>
<div class="note" id="o-n"></div>
<script>''' + COMMON + r'''
const RF=1e6,CB=10e-6;let sect=false,C={};const sectorOf=k=>"ACBD"[Math.floor(((k+12)%N)/24)];
function recompute(){const vab=+$('vab').value,vcd=+$('vcd').value,t=+$('tim').value;const rf=field({A:0,B:0,C:t,D:t});
const dc=Array.from({length:N},(_,k)=>sect?((sectorOf(k)==='A'||sectorOf(k)==='B')?vab:vcd):(vab+vcd)/2);C={rf,dc};
$('vab-o').textContent=vab+' V';$('vcd-o').textContent=vcd+' V';$('tim-o').textContent=t+'°';
$('o-dc').textContent=[0,24,48,72].map(k=>dc[k].toFixed(1)).join(' / ')+' V';let pk=0;for(let k=0;k<N;k++)pk=Math.max(pk,Math.abs(dc[k])+env(rf,k));$('o-pk').textContent=pk.toFixed(1)+' V';
const dv=vab-vcd;$('o-i').textContent=((sect?0:dv/RF)*1e6).toFixed(1)+' µA';$('o-p').textContent=(sect?'0':(dv*dv/RF*1e3).toFixed(4))+' mW';
$('o-n').textContent=sect?('Each blocking capacitor holds '+Math.abs(dv)+' V and stores '+(0.5*CB*dv*dv*1e3).toFixed(2)+' mJ, charged once from the supplies; no steady DC current.'):'The whole ring sits at the average of the two supplies; the higher supply drives a steady current through the 1 MΩ feeds into the lower one.';
$('m-cont').classList.toggle('on',!sect);$('m-sec').classList.toggle('on',sect);}
['vab','vcd','tim'].forEach(id=>$(id).addEventListener('input',recompute));$('m-cont').onclick=()=>{sect=false;recompute();};$('m-sec').onclick=()=>{sect=true;recompute();};
recompute();controls();
loop(()=>{let VS=8;for(let k=0;k<N;k++)VS=Math.max(VS,Math.abs(C.dc[k])+env(C.rf,k));VS=Math.ceil(VS/5)*5;
frame(VS,[[k=>C.dc[k],'#BA7517',1.5],[k=>C.dc[k]+val(C.rf,k,ph),'#1D9E75',2.5]],[k=>C.dc[k]+env(C.rf,k),k=>C.dc[k]-env(C.rf,k)],sect?[45,135,225,315]:[]);});
</script></main></body></html>
'''
for n, t in files.items():
    open(outdir + '/' + n, 'w').write(t)
print('ok', list(files))
