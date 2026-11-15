/* Offline benchmark report. Pure calculation helpers also run in Node tests. */
"use strict";
const BenchmarkReport = (() => {
  const stages = {construct:"Construction", solve:"Solve", fit:"Fitting", evaluate_501:"Evaluate 501 pressures", export:"Export"};
  const colours = ["#087e8b", "#4268cc", "#d47b45", "#8f5dad", "#4c9271", "#bf5470", "#74733b", "#576f92"];
  const stageColours = ["#87c9c3", "#087e8b", "#4268cc", "#d6a64f", "#a781b9"];
  const median = values => { const a = [...values].sort((a,b)=>a-b), n=a.length; return n%2 ? a[(n-1)/2] : (a[n/2-1]+a[n/2])/2; };
  const stable = value => JSON.stringify(value, (_,v) => v && !Array.isArray(v) && typeof v === "object" ? Object.fromEntries(Object.entries(v).sort(([a],[b])=>a.localeCompare(b))) : v);
  function validate(report) {
    if (report.schema_version !== 1 || !report.cases || !Object.keys(report.cases).length) throw Error("Expected a nonempty schema-version-1 benchmark report.");
    for (const [id,c] of Object.entries(report.cases)) {
      if (!c.case || !c.parameters || !Array.isArray(c.samples) || c.samples.length < 2) throw Error(`${id}: first invocation and repeated samples are required.`);
      for (const sample of c.samples) for (const key of Object.keys(stages)) {
        if (!Number.isFinite(sample[key]) || sample[key] < 0) throw Error(`${id}: invalid ${key} duration.`);
      }
    }
    return report;
  }
  function stats(c, stage) {
    const values = c.samples.map(s => stage === "total" ? Object.keys(stages).reduce((sum,k)=>sum+s[k],0) : s[stage]);
    const warm = values.slice(1);
    return {first:values[0], median:median(warm), minimum:Math.min(...warm), maximum:Math.max(...warm), repeats:warm.length};
  }
  function comparisonReason(a,b,id) {
    const first=a.cases[id],second=b.cases[id];
    if(!first||!second)return 'Missing case';
    if(stable(first.parameters)!==stable(second.parameters))return 'Different case parameters';
    if(first.raw_endpoints&&second.raw_endpoints){
      for(const key of Object.keys(first.raw_endpoints)){
        const x=first.raw_endpoints[key],y=second.raw_endpoints[key];
        if(!y||x.length!==y.length||x.some((v,i)=>!Number.isFinite(v)||!Number.isFinite(y[i])||Math.abs(v-y[i])>1e-5+1e-4*Math.abs(v)))return `Different endpoint ${key}`;
      }
    }
    return null;
  }
  const compatible = (a,b,id) => comparisonReason(a,b,id)===null;
  const escape = s => String(s).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const seconds = x => x === null || !Number.isFinite(x) ? "—" : x >= 1 ? `${x.toFixed(3)} s` : x >= .001 ? `${(x*1000).toFixed(2)} ms` : `${(x*1e6).toFixed(1)} µs`;
  function fitErrorSummary(fit) {
    const absolute = fit.absolute_error.filter(Number.isFinite).map(Math.abs);
    const relative = fit.relative_error_percent.filter(Number.isFinite).map(Math.abs);
    return {maximumAbsolute:absolute.length ? Math.max(...absolute) : null,
            maximumRelative:relative.length ? Math.max(...relative) : null};
  }
  function init() {
    const $ = id => document.getElementById(id);
    let runs = JSON.parse($("benchmark-data").textContent);
    const svg = (width,height,body,label) => `<svg role="img" aria-label="${escape(label)}" viewBox="0 0 ${width} ${height}" style="min-width:${width}px">${body}</svg>`;
    const text = (x,y,value,extra="") => `<text x="${x}" y="${y}" ${extra}>${escape(value)}</text>`;
    const line = (x1,y1,x2,y2,cls="grid") => `<line class="${cls}" x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}"/>`;
    const title = value => `<title>${escape(value)}</title>`;
    const shortName = c => `${c.case.fluid} · ${c.case.family.replaceAll("_"," ")} · η ${c.case.efficiency}`;
    function entries() {
      const ids = [...new Set(runs.flatMap(r=>Object.keys(r.report.cases)))];
      return ids.filter(id => $("mode").value === "all" || runs.find(r=>r.report.cases[id]).report.cases[id].case.mode === $("mode").value);
    }
    function value(run,id,stage,stat,ratio=false) {
      const c=run.report.cases[id];
      if (!c) return null;
      const result = stats(c,stage);
      if (!ratio) return {v:result[stat], min:result.minimum, max:result.maximum, repeats:result.repeats};
      if (!compatible(runs[0].report,run.report,id)) return null;
      const base=stats(runs[0].report.cases[id],stage)[stat];
      return base>0 ? {v:result[stat]/base,min:result.minimum/base,max:result.maximum/base,repeats:result.repeats} : null;
    }
    function grouped(ids) {
      const stage=$("stage").value, stat=$("stat").value, ratio=$("view").value === "ratio";
      const width=Math.max(1050,ids.length*Math.max(94,36*runs.length)+100), height=425, left=75, top=30, bottom=300;
      const maximum=Math.max(0,...ids.flatMap(id=>runs.map(run=>{const v=value(run,id,stage,stat,ratio);return v ? Math.max(v.v,stat==="median"?v.max:0) : 0;})));
      const ceiling=(maximum || 1)*1.12, y=v=>bottom-v/ceiling*(bottom-top);
      let body="";
      for(let i=0;i<=5;i++){const tick=ceiling*i/5;body+=line(left,y(tick),width-25,y(tick))+text(left-12,y(tick)+4,ratio?`${tick.toFixed(2)}×`:seconds(tick),'text-anchor="end"');}
      if(ratio && ceiling>=1) body+=line(left,y(1),width-25,y(1),"baseline");
      const group=(width-left-25)/Math.max(ids.length,1), bw=Math.min(38,(group-20)/runs.length);
      ids.forEach((id,index)=>{
        const center=left+group*(index+.5);
        runs.forEach((run,j)=>{
          const v=value(run,id,stage,stat,ratio); if(!v) return;
          const x=center+(j-runs.length/2)*bw;
          const formatted=ratio?`${v.v.toFixed(3)}×`:seconds(v.v);
          body+=`<rect class="bar case-bar" data-case="${escape(id)}" data-run="${j}" role="button" tabindex="0" aria-label="${escape(`Inspect ${run.label}: ${id}`)}" x="${x+2}" y="${y(v.v)}" width="${Math.max(1,bw-4)}" height="${bottom-y(v.v)}" rx="3" fill="${colours[j%colours.length]}">${title(`${run.label}\n${id}\n${formatted}\n${v.repeats} repeated measurements${stat==="median"?`\nRange: ${ratio?v.min.toFixed(3)+"–"+v.max.toFixed(3)+"×":seconds(v.min)+" – "+seconds(v.max)}`:""}`)}</rect>`;
          if(stat==="median") {const mid=x+bw/2;body+=`<g stroke="#344d5f">${line(mid,y(v.min),mid,y(v.max),"")}${line(mid-4,y(v.min),mid+4,y(v.min),"")}${line(mid-4,y(v.max),mid+4,y(v.max),"")}</g>`;}
          if(runs.length<4) body+=text(x+bw/2,y(stat==="median"?v.max:v.v)-9,formatted,'text-anchor="middle" class="value"');
        });
        const c=runs.find(r=>r.report.cases[id]).report.cases[id];
        body+=text(center+18,bottom+24,shortName(c),`transform="rotate(-28 ${center+18} ${bottom+24})" text-anchor="end"`);
      });
      $("groupChart").innerHTML=svg(width,height,body,"Grouped computational cost by process and measured run").replace('role="img"','role="group"');
      $("chartNote").textContent=(stat==="median"?"Bars: median of repeated invocations. Whiskers: observed minimum–maximum, not confidence intervals. ":"Bars: first invocation per case in the process. ")+(ratio?"Ratios divide by the baseline statistic; below 1× is faster. Missing cases and changed parameters have no ratio bar.":"All times use the same linear scale within this chart. Change the stage to inspect smaller costs.");
    }
    function stacked() {
      const id=$("detail").value, stat=$("stat").value;
      const available=runs.filter(r=>r.report.cases[id]);
      const width=1050, left=260, right=140, height=55+available.length*64;
      const max=Math.max(1e-12,...available.map(r=>stats(r.report.cases[id],"total")[stat]));
      let body="";
      available.forEach((run,index)=>{
        const c=run.report.cases[id], total=stats(c,"total")[stat];
        const samples=stat==="first"?[c.samples[0]]:c.samples.slice(1);
        const means=Object.keys(stages).map(key=>samples.reduce((sum,s)=>sum+s[key],0)/samples.length);
        const sum=means.reduce((a,b)=>a+b,0);
        let x=left; const y=20+index*64;
        body+=text(left-15,y+22,run.label,'text-anchor="end"');
        means.forEach((mean,j)=>{
          const contribution=sum?total*mean/sum:0, w=contribution/max*(width-left-right);
          body+=`<rect class="bar" x="${x}" y="${y}" width="${w}" height="34" fill="${stageColours[j]}">${title(`${Object.values(stages)[j]}\n${seconds(contribution)} contribution to displayed total\n${sum?(100*mean/sum).toFixed(1):0}% of mean stage cost`)}</rect>`;x+=w;
        });
        body+=text(x+12,y+22,seconds(total),'class="value"');
      });
      $("stackChart").innerHTML=svg(width,height,body,"Stacked computational stages by measured run");
    }
    let inspectedCase=null;
    function plot(id, series, xLabel, yLabel, options={}) {
      const width=id==='errorPlot'?1100:560, height=id==='errorPlot'?280:330, left=id==='errorPlot'?110:76, right=22, top=30, bottom=height-58;
      const points=series.flatMap(s=>s.x.map((x,i)=>[x,s.y[i]])).filter(([x,y])=>Number.isFinite(x)&&Number.isFinite(y));
      if(!points.length){$(id).innerHTML='<p class="notice">No finite samples available for this plot.</p>';return;}
      let xs=options.xRange||[Math.min(...points.map(p=>p[0])),Math.max(...points.map(p=>p[0]))];
      let ys=options.yRange||[Math.min(...points.map(p=>p[1])),Math.max(...points.map(p=>p[1]))];
      if(options.zero) ys=[Math.min(0,ys[0]),Math.max(0,ys[1])];
      const extent=(range)=>{const pad=Math.max((range[1]-range[0])*.08,Math.max(Math.abs(range[0]),Math.abs(range[1]))*1e-6,1e-9);return [range[0]-pad,range[1]+pad];};
      xs=extent(xs);ys=extent(ys);
      const px=x=>left+(x-xs[0])/(xs[1]-xs[0])*(width-left-right), py=y=>bottom-(y-ys[0])/(ys[1]-ys[0])*(bottom-top);
      const fmt=v=>Math.abs(v)<1e-12?'0':Number(v.toPrecision(5)).toString();
      let body=`<defs><clipPath id="clip-${id}"><rect x="${left}" y="${top}" width="${width-left-right}" height="${bottom-top}"/></clipPath><marker id="arrow-${id}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#087e8b"/></marker></defs>`;
      for(let i=0;i<=4;i++){
        const x=xs[0]+(xs[1]-xs[0])*i/4,y=ys[0]+(ys[1]-ys[0])*i/4;
        body+=line(px(x),top,px(x),bottom)+text(px(x),bottom+20,fmt(x),'text-anchor="middle"');
        body+=line(left,py(y),width-right,py(y))+text(left-9,py(y)+4,fmt(y),'text-anchor="end"');
      }
      body+=text((left+width-right)/2,height-12,xLabel,'text-anchor="middle"')+text(17,(top+bottom)/2,yLabel,`text-anchor="middle" transform="rotate(-90 17 ${(top+bottom)/2})"`);
      body+=`<g clip-path="url(#clip-${id})">`;
      for(const boundary of options.breakpoints||[]) body+=line(px(boundary),top,px(boundary),bottom,'baseline');
      if(options.zero) body+=line(left,py(0),width-right,py(0),'baseline');
      series.forEach(s=>{
        let d='',connected=false;
        const valid=[];
        s.x.forEach((x,i)=>{const y=s.y[i];if(!Number.isFinite(x)||!Number.isFinite(y)){connected=false;return;}valid.push([x,y]);d+=`${connected?'L':'M'}${px(x)},${py(y)} `;connected=true;});
        if(s.points){valid.forEach(([x,y])=>{body+=`<circle cx="${px(x)}" cy="${py(y)}" r="2" fill="${s.colour}" opacity=".6">${title(`${s.name}: ${fmt(x)}, ${fmt(y)}`)}</circle>`;});}
        else body+=`<path d="${d}" fill="none" stroke="${s.colour}" stroke-width="${s.width||2}" ${s.arrow?`marker-end="url(#arrow-${id})"`:''}>${title(s.name)}</path>`;
        if(s.arrow&&valid.length){const start=valid[0],end=valid.at(-1);for(const [point,name] of [[start,'Inlet'],[end,'Outlet']]){body+=`<circle cx="${px(point[0])}" cy="${py(point[1])}" r="4" fill="${name==='Inlet'?'white':s.colour}" stroke="${s.colour}" stroke-width="2">${title(name)}</circle>`+text(px(point[0])+7,py(point[1])+(name==='Inlet'?-9:17),name);}}
      });
      body+='</g>';$(id).innerHTML=svg(width,height,body,`${yLabel} against ${xLabel}`);
    }
    function inspector() {
      if(!inspectedCase)return;
      const run=runs[Number($('inspectRun').value)], c=run?.report.cases[inspectedCase];
      if(!c)return;
      $('inspectorTitle').textContent=shortName(c);
      $('inspectorSubtitle').textContent=`${run.label} · ${c.case.mode} · jaxprop ${run.report.jaxprop?.source_version||'unknown'}`;
      const d=c.diagnostics;
      if(!d){
        $('diagnosticNotice').innerHTML='<p class="notice">This historical run has no stored trajectory or reference-state curves. Replay it with --diagnostics-only using the original jaxprop wheel, or record a new benchmark run. Timing bars remain available.</p>';
        for(const id of ['tsPlot','fitPlot','errorPlot','inspectorMetrics','errorBasis','diagnosticProvenance'])$(id).replaceChildren();
        $('inspectProperty').replaceChildren();return;
      }
      $('diagnosticNotice').replaceChildren();
      const previous=$('inspectProperty').value;
      $('inspectProperty').innerHTML=Object.keys(d.fits).map(key=>`<option value="${escape(key)}">${escape(key.replaceAll('_',' '))}</option>`).join('');
      if(d.fits[previous])$('inspectProperty').value=previous;
      const variable=$('inspectProperty').value, f=d.fits[variable], metrics=fitErrorSummary(f);
      $('inspectorMetrics').textContent=`${d.trajectory.pressure.length} trajectory states · ${f.curve_pressure.length} polynomial samples · Maximum |relative error|: ${metrics.maximumRelative===null?'undefined at zero reference':metrics.maximumRelative.toPrecision(4)+'%'} · Maximum |absolute error|: ${metrics.maximumAbsolute===null?'unavailable':metrics.maximumAbsolute.toPrecision(4)+' '+f.unit}`;
      const trajectory=d.trajectory,sat=d.saturation;
      const ts=[{name:'Saturated liquid',x:sat.entropy_liquid.map(x=>x===null?null:x/1000),y:sat.temperature,colour:'#aab8c1'},
        {name:'Saturated vapor',x:sat.entropy_vapor.map(x=>x===null?null:x/1000),y:sat.temperature,colour:'#aab8c1'},
        {name:'Process trajectory',x:trajectory.entropy.map(x=>x===null?null:x/1000),y:trajectory.temperature,colour:'#087e8b',width:3,arrow:true}];
      const tsOptions=$('inspectZoom').value==='process'?{xRange:[Math.min(...ts[2].x.filter(Number.isFinite)),Math.max(...ts[2].x.filter(Number.isFinite))],yRange:[Math.min(...trajectory.temperature),Math.max(...trajectory.temperature)]}:{};
      plot('tsPlot',ts,'Specific entropy [kJ/(kg·K)]','Temperature [K]',tsOptions);
      const p=f.pressure.map(x=>x/1e5), smoothP=f.curve_pressure.map(x=>x/1e5), breaks=d.breakpoints.slice(1,-1).map(x=>x/1e5);
      plot('fitPlot',[{name:'ODE model states',x:p,y:f.reference,colour:'#087e8b',points:true},{name:'Polynomial fit',x:smoothP,y:f.curve_fitted,colour:'#d47b45'}],'Pressure [bar]',`${variable.replaceAll('_',' ')} [${f.unit}]`,{breakpoints:breaks});
      const relative=$('inspectError').value==='relative';
      plot('errorPlot',[{name:'Fit minus ODE state',x:p,y:relative?f.relative_error_percent:f.absolute_error,colour:'#4268cc',points:true}],'Pressure [bar]',relative?'Relative fit error [%]':`Absolute fit error [${f.unit}]`,{zero:true,breakpoints:breaks});
      $('errorBasis').textContent=d.error_basis+' Relative error is omitted where the reference is zero or negligible; use absolute error for those points.';
      $('diagnosticProvenance').textContent=JSON.stringify(d.provenance,null,2);
    }
    function openInspector(id,runIndex) {
      inspectedCase=id;
      $('inspectRun').innerHTML=runs.map((run,i)=>run.report.cases[id]?`<option value="${i}">${escape(run.label)}</option>`:'').join('');
      $('inspectRun').value=String(runIndex);inspector();
      if(!$('caseInspector').open)$('caseInspector').showModal();
    }
    $('groupChart').onclick=event=>{const bar=event.target.closest('.case-bar');if(bar)openInspector(bar.dataset.case,Number(bar.dataset.run));};
    $('groupChart').onkeydown=event=>{const bar=event.target.closest('.case-bar');if(bar&&['Enter',' '].includes(event.key)){event.preventDefault();openInspector(bar.dataset.case,Number(bar.dataset.run));}};
    $('closeInspector').onclick=()=>{$('caseInspector').close();};
    ['inspectRun','inspectProperty','inspectZoom','inspectError'].forEach(id=>$(id).onchange=inspector);
    function warnings() {
      const notices=[], baseline=runs[0].report;
      for(const run of runs){
        const failures=Object.keys(run.report.failures||{});
        if(failures.length)notices.push(`${run.label}: ${failures.length} benchmark case(s) failed. No timing bars are shown for failed cases. See run provenance for recorded errors.`);
      }
      for(const run of runs.slice(1)) {
        const r=run.report;
        if(r.platform!==baseline.platform || r.processor!==baseline.processor) notices.push(`${run.label}: platform or processor differs from the baseline; timing ratios may reflect machine differences.`);
        const changed=Object.keys(r.cases).filter(id=>baseline.cases[id]&&stable(baseline.cases[id].parameters)!==stable(r.cases[id].parameters));
        if(changed.length) notices.push(`${run.label}: ${changed.length} case(s) have different parameters. Their baseline ratios are omitted.`);
        const numerical=Object.keys(r.cases).filter(id=>comparisonReason(baseline,r,id)?.startsWith('Different endpoint'));
        if(numerical.length)notices.push(`${run.label}: ${numerical.length} case(s) have endpoint differences beyond the numerical comparison budget. Their timing ratios are omitted; absolute costs remain visible.`);
        const missing=Object.keys(baseline.cases).filter(id=>!r.cases[id]);
        if(missing.length) notices.push(`${run.label}: ${missing.length} baseline case(s) are missing. No zero-time bars are substituted.`);
        const deps=Object.keys(baseline.dependencies||{}).filter(k=>k!=="jaxprop"&&r.dependencies?.[k]!==baseline.dependencies[k]);
        if(deps.length) notices.push(`${run.label}: other dependency versions differ (${deps.join(", ")}); this is not an isolated jaxprop comparison.`);
      }
      $("warnings").innerHTML=notices.map(n=>`<div class="warning">${escape(n)}</div>`).join("");
    }
    function render() {
      runs.forEach(r=>validate(r.report));
      const ids=entries();
      $("runCount").textContent=runs.length;
      $("caseCount").textContent=ids.length;
      $("coverage").textContent=`${new Set(ids.map(id=>runs.find(r=>r.report.cases[id]).report.cases[id].case.fluid)).size} fluids · ${$("mode").value==="all"?"all physical models":$("mode").value}`;
      $("suiteCost").textContent=seconds(Object.values(runs[0].report.cases).reduce((sum,c)=>sum+stats(c,"total").median,0));
      $("singleNotice").hidden=runs.length>1;
      $("legend").replaceChildren();
      runs.forEach((run,i)=>{
        const key=document.createElement("div");key.className="key";
        const swatch=document.createElement("span");swatch.className="swatch";swatch.style.background=colours[i%colours.length];key.append(swatch);
        const input=document.createElement("input");input.value=run.label;input.setAttribute("aria-label",`Run ${i+1} label`);input.onchange=()=>{run.label=input.value.trim()||run.label;render();};key.append(input);
        if(i===0){const tag=document.createElement("span");tag.textContent="baseline";key.append(tag);}
        if(runs.length>1){const remove=document.createElement("button");remove.textContent="Remove";remove.onclick=()=>{runs.splice(i,1);render();};key.append(remove);}
        $("legend").append(key);
      });
      const detail=$("detail").value;
      $("detail").innerHTML=ids.map(id=>`<option value="${escape(id)}">${escape(id)}</option>`).join("");
      if(ids.includes(detail)) $("detail").value=detail;
      $("stageLegend").innerHTML=Object.values(stages).map((name,i)=>`<span class="key"><span class="swatch" style="background:${stageColours[i]}"></span>${escape(name)}</span>`).join("");
      grouped(ids);stacked();warnings();
      const stage=$("stage").value, stat=$("stat").value;
      $("ledgerNote").textContent=`${stage==="total"?"Total measured stages":stages[stage]} · ${stat==="median"?"repeated median":"first invocation"}. Ratios require matching case parameters. Ranges always describe repeated invocations.`;
      $("ledger").innerHTML=ids.flatMap(id=>runs.map((run,j)=>{
        const c=run.report.cases[id];if(!c)return "";
        const s=stats(c,stage), ratio=value(run,id,stage,stat,true);
        return `<tr><td><span class="swatch" style="background:${colours[j%colours.length]}"></span> ${escape(shortName(c))}<small>${escape(run.label)} · ${escape(c.case.mode)} · n=${s.repeats}</small></td><td>${seconds(s[stat])}</td><td>${seconds(s.minimum)} – ${seconds(s.maximum)}</td><td>${ratio?ratio.v.toFixed(3)+"×":"—"}</td><td>${escape((c.nfev||[]).join(" + "))}</td><td>${Number.isFinite(c.max_density_fit_relative_error)?(100*c.max_density_fit_relative_error).toFixed(4)+"%":"—"}</td></tr>`;
      })).join("");
      $("provenance").replaceChildren();
      runs.forEach(run=>{
        const detail=document.createElement("details"), summary=document.createElement("summary"), pre=document.createElement("pre");
        summary.textContent=`${run.label} · jaxprop ${run.report.jaxprop?.source_version||"unknown"} · ${run.filename}`;
        pre.className="meta";pre.textContent=JSON.stringify(Object.fromEntries(Object.entries(run.report).filter(([key])=>key!=="cases")),null,2);
        detail.append(summary,pre);$("provenance").append(detail);
      });
    }
    ["stage","stat","view","mode"].forEach(id=>$(id).onchange=render);
    $("detail").onchange=stacked;
    $("files").onchange=async event=>{
      $("error").textContent="";
      try {
        const pending=[];
        for(const file of event.target.files) {
          const report=validate(JSON.parse(await file.text()));
          let label=file.name.replace(/\.json$/i,"");
          const version=report.jaxprop?.source_version;if(version)label+=` · jaxprop ${version}`;
          pending.push({label,filename:file.name,report});
        }
        runs.push(...pending);render();
      } catch(error) { $("error").textContent=error.message; }
      event.target.value="";
    };
    $("save").onclick=()=>{
      const clone=document.documentElement.cloneNode(true);
      clone.querySelector('#caseInspector').removeAttribute('open');
      clone.querySelector("#benchmark-data").textContent=JSON.stringify(runs).replace(/</g,"\\u003c").replace(/>/g,"\\u003e").replace(/&/g,"\\u0026");
      const link=document.createElement("a"), url=URL.createObjectURL(new Blob(["<!doctype html>\n"+clone.outerHTML],{type:"text/html;charset=utf-8"}));
      link.href=url;link.download="barotropy-benchmark-report.html";link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
    };
    try { render(); } catch(error) { $("error").textContent=error.message; }
  }
  return {stages,median,stats,validate,compatible,comparisonReason,fitErrorSummary,init};
})();
if(typeof module!=="undefined") module.exports=BenchmarkReport;
if(typeof document!=="undefined") BenchmarkReport.init();
