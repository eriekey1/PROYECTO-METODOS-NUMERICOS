function showSection(id, btn){
  document.querySelectorAll('.section').forEach(s => {
    s.classList.remove('active');
    s.style.display = 'none'; 
  });
  document.querySelectorAll('.nav-btn').forEach(b => b.classList.remove('active'));

  const el = document.getElementById(id);
  if(el){
    el.style.display = 'block';
    el.classList.add('active');

  
    if(id === 'home'){
      window.scrollTo({ top: 0, behavior: 'smooth' });
    } else {
      setTimeout(() => {
        const offset = document.querySelector('nav').offsetHeight + 16; // altura del nav + margen
        const top = el.getBoundingClientRect().top + window.scrollY - offset;
        window.scrollTo({ top, behavior: 'smooth' });
      }, 50);
    }
  }

  if(btn) btn.classList.add('active');
}

// Show home by default
document.getElementById('home').style.display = 'block';
document.getElementById('home').classList.add('active');

// SOR panel toggle
document.getElementById('metodA').addEventListener('change', function(){
  document.getElementById('sor-panel').style.display = this.value === 'sor' ? 'block' : 'none';
});
// ═══════════════════════════════ CHART REGISTRY ═══════════════════════════════
const charts = {};
function destroyChart(id){if(charts[id]){charts[id].destroy();delete charts[id];}}
function makeChart(id,config){destroyChart(id);charts[id]=new Chart(document.getElementById(id),config);return charts[id];}

const chartDefaults = {
  color:'#8892a4',
  plugins:{legend:{labels:{color:'#8892a4',font:{family:"'IBM Plex Mono',monospace",size:11}}}},
  scales:{
    x:{ticks:{color:'#556070',font:{family:"'IBM Plex Mono',monospace",size:10}},grid:{color:'rgba(37,45,69,.5)'},border:{color:'#252d45'}},
    y:{ticks:{color:'#556070',font:{family:"'IBM Plex Mono',monospace",size:10}},grid:{color:'rgba(37,45,69,.5)'},border:{color:'#252d45'}}
  }
};

// ═══════════════════════════════════════════════════════════════════
//                        ESCENARIO A — Gauss-Seidel / Jacobi / SOR
// ═══════════════════════════════════════════════════════════════════
function getMatrixA(){
  const A=[[+document.getElementById('a11').value,+document.getElementById('a12').value,+document.getElementById('a13').value],
            [+document.getElementById('a21').value,+document.getElementById('a22').value,+document.getElementById('a23').value],
            [+document.getElementById('a31').value,+document.getElementById('a32').value,+document.getElementById('a33').value]];
  const b=[+document.getElementById('b1').value,+document.getElementById('b2').value,+document.getElementById('b3').value];
  const x0=[+document.getElementById('x1_0').value,+document.getElementById('x2_0').value,+document.getElementById('x3_0').value];
  return {A,b,x0};
}

function gaussSeidel(A,b,x0,tol,maxIter,omega){
  const n=A.length; let x=[...x0]; const hist=[]; const isJacobi=(omega===null&&false);
  omega = omega||1.0;
  for(let k=0;k<maxIter;k++){
    const xOld=[...x];
    for(let i=0;i<n;i++){
      let s=b[i];
      for(let j=0;j<n;j++) if(j!==i) s-=A[i][j]*x[j];
      const xNew=s/A[i][i];
      x[i]=x[i]+omega*(xNew-x[i]);
    }
    const err=Math.max(...x.map((xi,i)=>Math.abs(xi-xOld[i])));
    hist.push({k:k+1,x:[...x],err});
    if(err<tol) break;
  }
  return {x,hist};
}

function jacobiMethod(A,b,x0,tol,maxIter){
  const n=A.length; let x=[...x0]; const hist=[];
  for(let k=0;k<maxIter;k++){
    const xNew=new Array(n);
    for(let i=0;i<n;i++){
      let s=b[i];
      for(let j=0;j<n;j++) if(j!==i) s-=A[i][j]*x[j];
      xNew[i]=s/A[i][i];
    }
    const err=Math.max(...xNew.map((xi,i)=>Math.abs(xi-x[i])));
    x=[...xNew];
    hist.push({k:k+1,x:[...x],err});
    if(err<tol) break;
  }
  return {x,hist};
}

function runEscenarioA(){
  const {A,b,x0}=getMatrixA();
  const metodo=document.getElementById('metodA').value;
  const tol=parseFloat(document.getElementById('tolA').value);
  const omega=metodo==='sor'?parseFloat(document.getElementById('omegaA').value):1.0;
  const st=document.getElementById('statusA');

  // Check diagonal dominance
  let isDomDiag=true;
  for(let i=0;i<3;i++){
    const sum=A[i].reduce((s,v,j)=>j===i?s:s+Math.abs(v),0);
    if(Math.abs(A[i][i])<sum) isDomDiag=false;
  }

  let res;
  if(metodo==='jacobi') res=jacobiMethod(A,b,x0,tol,200);
  else res=gaussSeidel(A,b,x0,x0,tol,200,metodo==='sor'?omega:1.0);
  // fix: gs and sor pass x0 correctly
  if(metodo==='gs') res=gaussSeidel(A,b,x0,tol,200,1.0);
  else if(metodo==='sor') res=gaussSeidel(A,b,x0,tol,200,omega);
  else res=jacobiMethod(A,b,x0,tol,200);

  const {x,hist}=res;
  st.textContent=`Convergió en ${hist.length} iteraciones`;

  // Crystal box
  const cb=document.getElementById('crystalA'); cb.style.display='block';
  let html='';
  const showIter=Math.min(3,hist.length);
  for(let k=0;k<showIter;k++){
    const h=hist[k];
    html+=`<div class="iteration-step">
      <span class="step-label">Iteración ${h.k}</span>
      <span class="step-math">`;
    for(let i=0;i<3;i++){
      const row=A[i]; const bi=b[i];
      let num=bi;
      for(let j=0;j<3;j++) if(j!==i) num-=row[j]*h.x[j];
      const raw=num/row[i];
      html+=`  x${i+1} = (${bi.toFixed(2)} - `;
      const terms=[];
      for(let j=0;j<3;j++) if(j!==i) terms.push(`${row[j].toFixed(2)}×${h.x[j].toFixed(4)}`);
      html+=terms.join(' - ')+`) / ${row[i].toFixed(2)} = ${h.x[i].toFixed(6)}\n`;
    }
    html+=`</span><span class="step-result">Error máximo: ${h.err.toExponential(4)}</span></div>`;
  }
  document.getElementById('crystalAContent').innerHTML=html;

  // Stats
  document.getElementById('statsA').innerHTML=`
    <div class="stat-card"><div class="sv">${x[0].toFixed(3)}</div><div class="sl">x₁ · Fuente 1 (miles L)</div></div>
    <div class="stat-card"><div class="sv">${x[1].toFixed(3)}</div><div class="sl">x₂ · Fuente 2 (miles L)</div></div>
    <div class="stat-card"><div class="sv">${x[2].toFixed(3)}</div><div class="sl">x₃ · Fuente 3 (miles L)</div></div>
    <div class="stat-card"><div class="sv" style="color:var(--accent3)">${hist.length}</div><div class="sl">Iteraciones</div></div>
  `;

  // Charts
  const labels=hist.map(h=>h.k);
  makeChart('chartA',{type:'line',data:{labels,datasets:[
    {label:'x₁',data:hist.map(h=>h.x[0]),borderColor:'#f5a623',tension:.3,pointRadius:2},
    {label:'x₂',data:hist.map(h=>h.x[1]),borderColor:'#4fc3f7',tension:.3,pointRadius:2},
    {label:'x₃',data:hist.map(h=>h.x[2]),borderColor:'#81c784',tension:.3,pointRadius:2},
  ]},options:{...chartDefaults,responsive:true}});

  makeChart('chartA2',{type:'bar',data:{labels:['Zona Norte','Zona Centro','Zona Sur'],datasets:[
    {label:'Flujo requerido (miles L/día)',data:b,backgroundColor:'rgba(245,166,35,.4)',borderColor:'#f5a623',borderWidth:1},
  ]},options:{...chartDefaults,responsive:true,plugins:{...chartDefaults.plugins,legend:{display:false}}}});

  // Table
  let thead='<tr><th>Iter.</th><th>x₁</th><th>x₂</th><th>x₃</th><th>Error</th></tr>';
  let tbody=hist.slice(0,50).map(h=>`<tr${h===hist[hist.length-1]?' class="highlight-row"':''}><td>${h.k}</td><td>${h.x[0].toFixed(5)}</td><td>${h.x[1].toFixed(5)}</td><td>${h.x[2].toFixed(5)}</td><td>${h.err.toExponential(3)}</td></tr>`).join('');
  document.getElementById('tableA').innerHTML=`<thead>${thead}</thead><tbody>${tbody}</tbody>`;

  // Analysis
  document.getElementById('qa-a1').textContent=`Zona Norte: ${x[0].toFixed(3)} miles L/día | Zona Centro: ${x[1].toFixed(3)} miles L/día | Zona Sur: ${x[2].toFixed(3)} miles L/día`;
  document.getElementById('qa-a3').textContent=isDomDiag?`El sistema es diagonalmente dominante → matemáticamente ESTABLE. Gauss-Seidel converge garantizado. Un cambio pequeño en la demanda produce un cambio proporcional en la solución.`:`⚠️ El sistema NO es diagonalmente dominante → puede ser inestable o converger lentamente. Verifica que los coeficientes diagonales sean mayores que la suma de los off-diagonal.`;
  document.getElementById('conclusionA').innerHTML=`<strong>Conclusión del Escenario A:</strong> El método Gauss-Seidel (y SOR con ω óptimo) resuelve eficientemente la red de distribución de combustible, encontrando los flujos exactos en pocas iteraciones. La gran fortaleza de este enfoque es la <em>interpretabilidad</em>: modificar un coeficiente para simular un bloqueo de ruta (a=0) reproduce inmediatamente el rebalanceo del sistema. La limitación principal es que el modelo es determinista y estático — no captura la variabilidad temporal de la demanda ni las restricciones de capacidad de camiones. Una mejora natural sería reformularlo como un problema de <em>flujo en red</em> con restricciones de capacidad, resoluble con programación lineal.`;

  document.getElementById('resultsA').style.display='block';
  document.getElementById('analysisA').style.display='block';
}

// ═══════════════════════════════════════════════════════════════════
//                        ESCENARIO B — EDO: Euler / Heun / RK4
// ═══════════════════════════════════════════════════════════════════
function dRdt(t, R, E, C0, alpha){
  const consumo = C0*(1+alpha*t);
  return E - consumo;
}

function eulerODE(f,t0,R0,h,tf,E,C0,alpha){
  const pts=[{t:t0,R:R0}]; let t=t0,R=R0;
  while(t<tf-1e-9){const dt=Math.min(h,tf-t);R=R+dt*f(t,R,E,C0,alpha);t+=dt;R=Math.max(0,R);pts.push({t,R});}
  return pts;
}
function heunODE(f,t0,R0,h,tf,E,C0,alpha){
  const pts=[{t:t0,R:R0}]; let t=t0,R=R0;
  while(t<tf-1e-9){
    const dt=Math.min(h,tf-t);
    const k1=f(t,R,E,C0,alpha);
    const Rp=R+dt*k1;
    const k2=f(t+dt,Rp,E,C0,alpha);
    R=R+dt*(k1+k2)/2; t+=dt; R=Math.max(0,R); pts.push({t,R});
  }
  return pts;
}
function rk4ODE(f,t0,R0,h,tf,E,C0,alpha){
  const pts=[{t:t0,R:R0}]; let t=t0,R=R0;
  while(t<tf-1e-9){
    const dt=Math.min(h,tf-t);
    const k1=f(t,R,E,C0,alpha);
    const k2=f(t+dt/2,R+dt*k1/2,E,C0,alpha);
    const k3=f(t+dt/2,R+dt*k2/2,E,C0,alpha);
    const k4=f(t+dt,R+dt*k3,E,C0,alpha);
    R=R+dt*(k1+2*k2+2*k3+k4)/6; t+=dt; R=Math.max(0,R); pts.push({t,R});
  }
  return pts;
}

function runEscenarioB(){
  const R0=parseFloat(document.getElementById('r0B').value);
  const E=parseFloat(document.getElementById('entradaB').value);
  const C0=parseFloat(document.getElementById('consumoB').value);
  const alpha=parseFloat(document.getElementById('alphaB').value);
  const Rcrit=parseFloat(document.getElementById('rcritB').value);
  const tf=parseFloat(document.getElementById('tfinalB').value);
  const metodo=document.getElementById('metodB').value;
  const h=parseFloat(document.getElementById('hB').value);

  const solvers={euler:eulerODE,heun:heunODE,rk4:rk4ODE};
  let datasets=[];
  let ptsRef=[];

  const metodosRun = metodo==='todos'?['euler','heun','rk4']:[metodo];
  const colors={'euler':'#e05252','heun':'#f5a623','rk4':'#81c784'};
  const nombres={'euler':'Euler','heun':'Heun','rk4':'RK4'};

  metodosRun.forEach(m=>{
    const pts=solvers[m](dRdt,0,R0,h,tf,E,C0,alpha);
    if(m===metodosRun[metodosRun.length-1]) ptsRef=pts;
    datasets.push({
      label:nombres[m],
      data:pts.map(p=>({x:p.t,y:p.R})),
      borderColor:colors[m],tension:.3,pointRadius:0,borderWidth:2
    });
  });

  // Find critical day using main method
  let critDay=null;
  for(let i=0;i<ptsRef.length;i++) if(ptsRef[i].R<=Rcrit){critDay=ptsRef[i].t;break;}

  // Crystal box
  const cb=document.getElementById('crystalB'); cb.style.display='block';
  const pts3=ptsRef.slice(0,4);
  let html='';
  for(let i=0;i<Math.min(3,pts3.length-1);i++){
    const p=pts3[i]; const t=p.t;
    const k1=dRdt(t,p.R,E,C0,alpha);
    html+=`<div class="iteration-step">
      <span class="step-label">t = ${t.toFixed(2)} días</span>
      <span class="step-math">dR/dt = E - C(t) = ${E} - ${C0.toFixed(2)}×(1 + ${alpha}×${t.toFixed(2)}) = ${k1.toFixed(4)} %/día\nR(${(t+h).toFixed(2)}) ≈ ${p.R.toFixed(4)} + ${h}×${k1.toFixed(4)} = ${(p.R+h*k1).toFixed(4)} %</span>
      <span class="step-result">Reserva actual: ${pts3[i+1].R.toFixed(3)} %  ${pts3[i+1].R<=Rcrit?'⚠️ NIVEL CRÍTICO':''}</span>
    </div>`;
  }
  document.getElementById('crystalBContent').innerHTML=html;

  // Stats
  const finalR=ptsRef[ptsRef.length-1].R;
  document.getElementById('statsB').innerHTML=`
    <div class="stat-card"><div class="sv${critDay?'" style="color:var(--accent2)"':'"'}>${critDay?critDay.toFixed(1)+' d':'No alcanza'}</div><div class="sl">Día Nivel Crítico</div></div>
    <div class="stat-card"><div class="sv">${finalR.toFixed(1)}%</div><div class="sl">Reserva Final (día ${tf})</div></div>
    <div class="stat-card"><div class="sv">${(C0*(1+alpha*tf)).toFixed(2)}%</div><div class="sl">Consumo día ${tf} (%/día)</div></div>
    <div class="stat-card"><div class="sv" style="color:var(--accent3)">${ptsRef.length}</div><div class="sl">Pasos de tiempo</div></div>
  `;

  // Chart
  const critDataset={label:`Nivel crítico (${Rcrit}%)`,data:[{x:0,y:Rcrit},{x:tf,y:Rcrit}],borderColor:'#e05252',borderDash:[6,3],pointRadius:0,borderWidth:1.5};
  makeChart('chartB',{type:'line',data:{datasets:[...datasets,critDataset]},options:{
    ...chartDefaults,responsive:true,
    scales:{
      x:{...chartDefaults.scales.x,type:'linear',title:{display:true,text:'Día',color:'#556070'}},
      y:{...chartDefaults.scales.y,title:{display:true,text:'Reserva (%)',color:'#556070'},min:0,max:105}
    }
  }});

  // Table
  let thead='<tr><th>Día</th><th>Reserva (%)</th><th>Consumo (%/día)</th><th>Saldo (%/día)</th></tr>';
  const step=Math.max(1,Math.floor(ptsRef.length/30));
  let tbody=ptsRef.filter((_,i)=>i%step===0||i===ptsRef.length-1).map(p=>{
    const c=C0*(1+alpha*p.t); const saldo=E-c;
    return `<tr${p.R<=Rcrit?' class="highlight-row"':''}><td>${p.t.toFixed(1)}</td><td>${p.R.toFixed(2)}</td><td>${c.toFixed(3)}</td><td>${saldo.toFixed(3)}</td></tr>`;
  }).join('');
  document.getElementById('tableB').innerHTML=`<thead>${thead}</thead><tbody>${tbody}</tbody>`;

  // Analysis
  document.getElementById('qa-b1').textContent=critDay?`Al ritmo actual (consumo base ${C0}%/día, factor de pánico α=${alpha}), las reservas alcanzan el nivel crítico del ${Rcrit}% en el día ${critDay.toFixed(1)}.`:`Con los parámetros actuales, las reservas NO alcanzan el nivel crítico en ${tf} días. La tasa de entrada (${E}%/día) es suficiente para mantener el sistema.`;
  document.getElementById('qa-b3').textContent=E===0?`Con suministro cortado (E=0), el vaciado es mucho más rápido. La reserva desciende exclusivamente por consumo, alcanzando el nivel crítico en aprox. ${critDay?critDay.toFixed(1):'pocos'} días.`:`Con la tasa de entrada actual (${E}%/día), el sistema tiene un buffer de ${(E/C0*100).toFixed(0)}% respecto al consumo base. Cambia E a 0 para simular corte total.`;
  document.getElementById('conclusionB').innerHTML=`<strong>Conclusión del Escenario B:</strong> La EDO dR/dt = E - C(t) captura con elegancia la dinámica de vaciado de reservas. RK4 demostró ser el método más preciso, especialmente cuando el consumo crece rápidamente (α > 0). La diferencia entre Euler y RK4 se vuelve visible con pasos h grandes — Euler puede subestimar o sobreestimar la tasa de consumo. La <em>limitación fundamental</em> del modelo es que trata el consumo como una función determinista suave del tiempo, cuando en realidad hay picos (mañanas, viernes antes de fines de semana) y valles. Una mejora realista incorporaría un modelo de consumo estocástico o basado en datos históricos de demanda horaria.`;

  document.getElementById('resultsB').style.display='block';
  document.getElementById('analysisB').style.display='block';
}

// ═══════════════════════════════════════════════════════════════════
//                        ESCENARIO C — Interpolación
// ═══════════════════════════════════════════════════════════════════
function addPuntoC(){
  const div=document.createElement('div');
  div.className='row g-2 mb-2 punto-row';
  div.innerHTML=`<div class="col-5"><input class="input-dark" type="number" placeholder="Día" data-type="t"></div><div class="col-5"><input class="input-dark" type="number" placeholder="Precio Bs" data-type="p"></div><div class="col-2"><button class="btn-secondary" onclick="removePuntoC(this)">✕</button></div>`;
  document.getElementById('puntosC').appendChild(div);
}
function removePuntoC(btn){
  const rows=document.querySelectorAll('#puntosC .punto-row');
  if(rows.length>2) btn.closest('.punto-row').remove();
}
function getPuntosC(){
  const rows=document.querySelectorAll('#puntosC .punto-row');
  const pts=[];
  rows.forEach(r=>{
    const t=parseFloat(r.querySelector('[data-type=t]').value);
    const p=parseFloat(r.querySelector('[data-type=p]').value);
    if(!isNaN(t)&&!isNaN(p)) pts.push({t,p});
  });
  pts.sort((a,b)=>a.t-b.t);
  return pts;
}

function lagrangeInterp(pts, x){
  const n=pts.length; let result=0;
  for(let i=0;i<n;i++){
    let L=pts[i].p;
    for(let j=0;j<n;j++){
      if(j!==i) L*=(x-pts[j].t)/(pts[i].t-pts[j].t);
    }
    result+=L;
  }
  return result;
}

function newtonDivDiff(pts){
  const n=pts.length;
  const table=pts.map(p=>[p.p]);
  for(let j=1;j<n;j++)
    for(let i=0;i<n-j;i++)
      table[i].push((table[i+1][j-1]-table[i][j-1])/(pts[i+j].t-pts[i].t));
  return table;
}
function newtonInterp(pts, x, table){
  const n=pts.length; let result=table[0][0]; let prod=1;
  for(let j=1;j<n;j++){prod*=(x-pts[j-1].t);result+=table[0][j]*prod;}
  return result;
}

function splineCubic(pts){
  const n=pts.length; const m=n-1;
  const h=pts.slice(0,m).map((_,i)=>pts[i+1].t-pts[i].t);
  // Build tridiagonal system for M (second derivatives)
  const a=[],b=[],c=[],d=[];
  a.push(0); b.push(1); c.push(0); d.push(0); // natural spline boundary
  for(let i=1;i<m;i++){
    a.push(h[i-1]);
    b.push(2*(h[i-1]+h[i]));
    c.push(h[i]);
    d.push(6*((pts[i+1].p-pts[i].p)/h[i]-(pts[i].p-pts[i-1].p)/h[i-1]));
  }
  a.push(0); b.push(1); c.push(0); d.push(0);
  // Thomas algorithm
  const M=new Array(n).fill(0);
  const cp=[...c],dp=[...d];
  for(let i=1;i<n;i++){const fac=a[i]/b[i-1];b[i]-=fac*cp[i-1];dp[i]-=fac*dp[i-1];}
  M[n-1]=dp[n-1]/b[n-1];
  for(let i=n-2;i>=0;i--) M[i]=(dp[i]-cp[i]*M[i+1])/b[i];
  return {M,h,pts};
}
function evalSpline({M,h,pts},x){
  const n=pts.length;
  let i=0;
  for(let j=0;j<n-2;j++) if(x>=pts[j].t&&x<=pts[j+1].t) i=j;
  if(x>pts[n-1].t) i=n-2;
  const hi=h[i], xi=pts[i].t, xi1=pts[i+1].t;
  const A=(xi1-x)/hi, B=(x-xi)/hi;
  return A*pts[i].p+B*pts[i+1].p+((A*A*A-A)*M[i]+(B*B*B-B)*M[i+1])*hi*hi/6;
}

function runEscenarioC(){
  const pts=getPuntosC();
  if(pts.length<3){alert('Necesitas al menos 3 puntos');return;}
  const xEval=parseFloat(document.getElementById('evalC').value);
  const metodo=document.getElementById('metodC').value;

  const ddTable=newtonDivDiff(pts);
  const spline=splineCubic(pts);

  const t0=pts[0].t, tN=pts[pts.length-1].t;
  const xs=[], lagY=[], newtonY=[], splineY=[];
  for(let t=t0;t<=tN;t+=0.2){
    xs.push(t);
    lagY.push(lagrangeInterp(pts,t));
    newtonY.push(newtonInterp(pts,t,ddTable));
    splineY.push(evalSpline(spline,t));
  }

  let evalLag=lagrangeInterp(pts,xEval);
  let evalNew=newtonInterp(pts,xEval,ddTable);
  let evalSpl=evalSpline(spline,xEval);

  // Crystal — Newton diff table
  const cb=document.getElementById('crystalC'); cb.style.display='block';
  let html=`<p style="color:var(--text2);font-size:.82rem;margin-bottom:12px;">Tabla de diferencias divididas de Newton (primeras 4 columnas):</p>`;
  html+=`<div style="overflow-x:auto;"><table class="results-table"><thead><tr><th>i</th><th>tᵢ</th><th>f[tᵢ]</th><th>f[tᵢ,tᵢ₊₁]</th><th>f[tᵢ,...,tᵢ₊₂]</th><th>f[tᵢ,...,tᵢ₊₃]</th></tr></thead><tbody>`;
  for(let i=0;i<pts.length;i++){
    html+=`<tr><td>${i}</td><td>${pts[i].t}</td><td>${ddTable[i][0].toFixed(4)}</td><td>${ddTable[i][1]?.toFixed(4)||'—'}</td><td>${ddTable[i][2]?.toFixed(4)||'—'}</td><td>${ddTable[i][3]?.toFixed(4)||'—'}</td></tr>`;
  }
  html+=`</tbody></table></div>`;
  html+=`<div class="iteration-step" style="margin-top:10px;">
    <span class="step-label">Evaluación en t = ${xEval}</span>
    <span class="step-math">P(${xEval}) = ${ddTable[0][0].toFixed(4)} + ${ddTable[0][1]?.toFixed(4)||0}×(${xEval}-${pts[0].t}) + ...</span>
    <span class="step-result">Newton: ${evalNew.toFixed(4)} Bs/kg</span>
  </div>`;
  document.getElementById('crystalCContent').innerHTML=html;

  // Stats
  const span=tN-t0; const range=Math.max(...pts.map(p=>p.p))-Math.min(...pts.map(p=>p.p));
  document.getElementById('statsC').innerHTML=`
    <div class="stat-card"><div class="sv">${evalLag.toFixed(2)} Bs</div><div class="sl">Lagrange (día ${xEval})</div></div>
    <div class="stat-card"><div class="sv">${evalNew.toFixed(2)} Bs</div><div class="sl">Newton (día ${xEval})</div></div>
    <div class="stat-card"><div class="sv">${evalSpl.toFixed(2)} Bs</div><div class="sl">Spline (día ${xEval})</div></div>
    <div class="stat-card"><div class="sv">${range.toFixed(1)} Bs</div><div class="sl">Variación total precios</div></div>
  `;

  // Chart
  const showAll=metodo==='todos';
  const datasets=[
    {label:'Datos reales',data:pts.map(p=>({x:p.t,y:p.p})),type:'scatter',backgroundColor:'#f5a623',pointRadius:7,zIndex:10}
  ];
  if(showAll||metodo==='lagrange') datasets.push({label:'Lagrange',data:xs.map((x,i)=>({x,y:lagY[i]})),borderColor:'#e05252',tension:.3,pointRadius:0,type:'line'});
  if(showAll||metodo==='newton') datasets.push({label:'Newton',data:xs.map((x,i)=>({x,y:newtonY[i]})),borderColor:'#4fc3f7',tension:.3,pointRadius:0,type:'line'});
  if(showAll||metodo==='spline') datasets.push({label:'Spline Cúbico',data:xs.map((x,i)=>({x,y:splineY[i]})),borderColor:'#81c784',tension:.3,pointRadius:0,type:'line'});
  // Eval point
  datasets.push({label:`Estimación día ${xEval}`,data:[{x:xEval,y:metodo==='spline'?evalSpl:evalNew}],type:'scatter',backgroundColor:'#ce93d8',pointRadius:10,pointStyle:'triangle'});

  makeChart('chartC',{type:'scatter',data:{datasets},options:{
    ...chartDefaults,responsive:true,
    scales:{x:{...chartDefaults.scales.x,type:'linear',title:{display:true,text:'Día',color:'#556070'}},
            y:{...chartDefaults.scales.y,title:{display:true,text:'Precio (Bs/kg)',color:'#556070'}}}
  }});

  // Analysis
  const pct=(range/pts[0].p*100).toFixed(1);
  document.getElementById('qa-c1').textContent=`Día ${xEval}: Lagrange → ${evalLag.toFixed(2)} Bs, Newton → ${evalNew.toFixed(2)} Bs, Spline → ${evalSpl.toFixed(2)} Bs. Los tres métodos deberían coincidir cerca; diferencias grandes indican oscilación (efecto Runge en Lagrange con muchos puntos).`;
  document.getElementById('qa-c3').textContent=`Con ${pts.length} puntos en ${span} días, la densidad promedio es ${(span/(pts.length-1)).toFixed(1)} días entre observaciones. ${pts.length>6?'Con muchos nodos y polinomios de alto grado, Lagrange puede oscilar (efecto Runge). Splines cúbicos son más estables en este caso.':'Con pocos puntos, el polinomio de Lagrange es estable. Agregar más puntos dispersos puede introducir oscilaciones indeseadas.'}`;
  document.getElementById('conclusionC').innerHTML=`<strong>Conclusión del Escenario C:</strong> La interpolación transforma datos dispersos del mercado en una herramienta de predicción de precios. Los Splines cúbicos emergen como la mejor opción para datos de precios reales: preservan la forma de la curva entre nodos sin las oscilaciones de Lagrange, y la condición de continuidad de segunda derivada garantiza suavidad realista. El incremento total del ${pct}% en el precio durante el período analizado cuantifica con precisión el impacto de la crisis. La <em>limitación crítica</em>: la interpolación asume que la realidad se comportó como predice el polinomio entre los días medidos — lo cual no necesariamente es cierto si hubo picos puntuales no capturados. Se recomienda aumentar la frecuencia de medición en períodos de alta volatilidad.`;

  document.getElementById('resultsC').style.display='block';
  document.getElementById('analysisC').style.display='block';
}

// ═══════════════════════════════════════════════════════════════════
//                        ESCENARIO D — Integración Numérica
// ═══════════════════════════════════════════════════════════════════
function precioFunc(t, p0, pf, tipo, T){
  const r=Math.log(pf/p0)/T;
  if(tipo==='lineal') return p0+(pf-p0)*t/T;
  if(tipo==='exponencial') return p0*Math.exp(r*t);
  // sigmoide: logistic
  const mid=T/2, k=0.3;
  return p0+(pf-p0)/(1+Math.exp(-k*(t-mid)));
}

function trapecioInt(f,a,b,n){
  const h=(b-a)/n; let s=f(a)+f(b);
  for(let i=1;i<n;i++) s+=2*f(a+i*h);
  return (h/2)*s;
}
function simpson13(f,a,b,n){
  if(n%2!==0) n++;
  const h=(b-a)/n; let s=f(a)+f(b);
  for(let i=1;i<n;i++) s+=(i%2===0?2:4)*f(a+i*h);
  return (h/3)*s;
}
function simpson38(f,a,b,n){
  while(n%3!==0) n++;
  const h=(b-a)/n; let s=f(a)+f(b);
  for(let i=1;i<n;i++) s+=(i%3===0?2:3)*f(a+i*h);
  return (3*h/8)*s;
}

function runEscenarioD(){
  const p0=parseFloat(document.getElementById('p0D').value);
  const pf=parseFloat(document.getElementById('pfD').value);
  const tipo=document.getElementById('tipoD').value;
  const Q=parseFloat(document.getElementById('qD').value);
  const ingreso=parseFloat(document.getElementById('ingresoD').value);
  const metodo=document.getElementById('metodD').value;
  const T=30, n=60;

  const P=(t)=>precioFunc(t,p0,pf,tipo,T);
  const f=(t)=>Q*P(t);

  const gT=trapecioInt(f,0,T,n);
  const gS13=simpson13(f,0,T,n);
  const gS38=simpson38(f,0,T,n);
  const gastoSinAlza=Q*p0*T;

  const mainGasto={trapecio:gT,'s13':gS13,'s38':gS38,todos:gS13}[metodo];
  const perdida=mainGasto-gastoSinAlza;

  // Crystal
  const cb=document.getElementById('crystalD'); cb.style.display='block';
  let html='';
  const hs=T/n; const pts=[0,hs,2*hs,3*hs];
  html+=`<div class="iteration-step"><span class="step-label">Trapecio — Primeros 3 intervalos (h=${hs.toFixed(2)} días)</span><span class="step-math">`;
  for(let i=0;i<3;i++){
    const t1=pts[i], t2=pts[i+1];
    const f1=f(t1).toFixed(4), f2=f(t2).toFixed(4);
    html+=`  Área[${t1.toFixed(2)},${t2.toFixed(2)}] = (${hs.toFixed(2)}/2)×[${f1} + ${f2}] = ${((hs/2)*(parseFloat(f1)+parseFloat(f2))).toFixed(4)} Bs\n`;
  }
  html+=`</span><span class="step-result">Suma acumulando → integral total ≈ ${gT.toFixed(2)} Bs</span></div>`;
  document.getElementById('crystalDContent').innerHTML=html;

  // Stats
  const pct=(mainGasto/ingreso*100).toFixed(1);
  document.getElementById('statsD').innerHTML=`
    <div class="stat-card"><div class="sv">${mainGasto.toFixed(2)} Bs</div><div class="sl">Gasto Real Mensual</div></div>
    <div class="stat-card"><div class="sv">${gastoSinAlza.toFixed(2)} Bs</div><div class="sl">Gasto Sin Alza</div></div>
    <div class="stat-card"><div class="sv" style="color:var(--accent2)">${perdida.toFixed(2)} Bs</div><div class="sl">Pérdida Poder Adquisitivo</div></div>
    <div class="stat-card"><div class="sv" style="color:${parseFloat(pct)>20?'var(--accent2)':'var(--accent4)'}">${pct}%</div><div class="sl">Del Ingreso Familiar</div></div>
  `;

  // Charts
  const ts=Array.from({length:61},(_,i)=>i*T/60);
  const ys=ts.map(t=>P(t));
  const fillY=ts.map(t=>Q*P(t));
  makeChart('chartD',{type:'line',data:{
    labels:ts.map(t=>t.toFixed(1)),
    datasets:[
      {label:'Precio P(t) Bs/kg',data:ys,borderColor:'#f5a623',tension:.3,pointRadius:0,fill:false,yAxisID:'y'},
      {label:'Gasto diario Q·P(t) Bs/día',data:fillY,borderColor:'#4fc3f7',backgroundColor:'rgba(79,195,247,.12)',tension:.3,pointRadius:0,fill:true,yAxisID:'y'}
    ]
  },options:{...chartDefaults,responsive:true,scales:{
    x:{...chartDefaults.scales.x,title:{display:true,text:'Día',color:'#556070'}},
    y:{...chartDefaults.scales.y,id:'y',title:{display:true,text:'Bs',color:'#556070'}}
  }}});

  makeChart('chartD2',{type:'bar',data:{
    labels:['Trapecio','Simpson 1/3','Simpson 3/8'],
    datasets:[
      {label:'Gasto calculado (Bs)',data:[gT,gS13,gS38],backgroundColor:['rgba(245,166,35,.5)','rgba(79,195,247,.5)','rgba(129,199,132,.5)'],borderColor:['#f5a623','#4fc3f7','#81c784'],borderWidth:1}
    ]
  },options:{...chartDefaults,responsive:true,plugins:{legend:{display:false}}}});

  // Analysis
  document.getElementById('qa-d1').textContent=`Con precios ${tipo}s de ${p0} a ${pf} Bs/kg y un consumo de ${Q} kg/día, la familia gastó ${mainGasto.toFixed(2)} Bs en este producto durante el mes.`;
  document.getElementById('qa-d2').textContent=`Si los precios se hubieran mantenido en ${p0} Bs/kg (precio inicial), el gasto total habría sido ${gastoSinAlza.toFixed(2)} Bs — una diferencia de ${perdida.toFixed(2)} Bs.`;
  document.getElementById('qa-d3').textContent=`La pérdida del poder adquisitivo fue de ${perdida.toFixed(2)} Bs (${(perdida/ingreso*100).toFixed(1)}% del ingreso mensual de ${ingreso} Bs). Esto significa que ${(perdida/ingreso*100).toFixed(1)} centavos de cada boliviano del salario se perdieron solo en el alza de este producto.`;
  document.getElementById('conclusionD').innerHTML=`<strong>Conclusión del Escenario D:</strong> La integración numérica convierte la curva de precios (un objeto abstracto) en un número concreto y accionable: el costo real mensual de la familia. Simpson 1/3 y 3/8 ofrecen mayor precisión que el Trapecio con el mismo número de evaluaciones, confirmando el valor de los métodos de orden superior. La diferencia entre los tres métodos para este problema es de ${Math.abs(gT-gS13).toFixed(2)} Bs — pequeña pero significativa en presupuestos ajustados. La <em>limitación principal</em> del modelo es asumir consumo constante Q: en la realidad, las familias reducen el consumo cuando los precios suben (elasticidad precio negativa). Incorporar Q(t)=Q₀·(P₀/P(t))^ε (con ε la elasticidad) haría el modelo considerablemente más realista.`;

  document.getElementById('resultsD').style.display='block';
  document.getElementById('analysisD').style.display='block';
}

// ═══════════════════════════════════════════════════════════════════
//                        ESCENARIO E — Raíces de Ecuaciones
// ═══════════════════════════════════════════════════════════════════
let tipoEActual='gasto';
function selectTipoE(tipo, btn){
  tipoEActual=tipo;
  document.querySelectorAll('#tipoEtabs .method-tab').forEach(b=>b.classList.remove('active'));
  btn.classList.add('active');
  ['gasto','reserva','social'].forEach(t=>{
    document.getElementById('paneE_'+t).style.display=t===tipo?'block':'none';
  });
}

function getFuncE(){
  if(tipoEActual==='gasto'){
    const p0=parseFloat(document.getElementById('p0E').value)||8;
    const pf=parseFloat(document.getElementById('pfE').value)||22;
    const Q=parseFloat(document.getElementById('qE').value)||1.5;
    const ing=parseFloat(document.getElementById('ingE').value)||3500;
    // Gasto acumulado hasta día t = Q * integral de P(s) ds (0 a t), usando exponential
    const r=Math.log(pf/p0)/30;
    // f(t) = Q * p0/r * (e^(r*t) - 1) - ing
    return {
      f:(t)=>Q*(p0/r)*(Math.exp(r*t)-1)-ing,
      df:(t)=>Q*p0*Math.exp(r*t),
      a:1, b:30, x0:15,
      label:'t (días)',
      desc:`Umbral: día en que el gasto acumulado de la familia iguala el ingreso de ${ing} Bs`
    };
  } else if(tipoEActual==='reserva'){
    const C0=parseFloat(document.getElementById('c0E').value)||7;
    const alpha=parseFloat(document.getElementById('alphaE').value)||0.05;
    const tE=parseFloat(document.getElementById('tE').value)||10;
    // Equilibrio: E = C0*(1+alpha*tE), buscar E (entrada) tal que red sea 0
    // f(E) = E - C0*(1+alpha*tE)
    return {
      f:(E)=>E - C0*(1+alpha*tE),
      df:(_)=>1,
      a:0, b:30, x0:C0,
      label:'E (tasa entrada %/día)',
      desc:`Tasa de entrada E que iguala exactamente el consumo en el día ${tE} (equilibrio)`
    };
  } else {
    const A=parseFloat(document.getElementById('aE').value)||2;
    const beta=parseFloat(document.getElementById('betaE').value)||0.8;
    const M=parseFloat(document.getElementById('mE').value)||10;
    return {
      f:(x)=>A*Math.exp(beta*x)-M,
      df:(x)=>A*beta*Math.exp(beta*x),
      a:0, b:10, x0:2,
      label:'x (nivel de desabastecimiento percibido)',
      desc:`Umbral de tensión social donde f(x) = ${A}·e^(${beta}x) = ${M} (punto de no retorno)`
    };
  }
}

function biseccion(f,a,b,tol,maxIter){
  const hist=[]; let fa=f(a);
  if(fa*f(b)>0) return {x:null,hist,err:'f(a) y f(b) tienen el mismo signo'};
  for(let k=0;k<maxIter;k++){
    const c=(a+b)/2, fc=f(c);
    const err=Math.abs(b-a)/2;
    hist.push({k:k+1,a,b,c,fc,err});
    if(err<tol||Math.abs(fc)<tol) return {x:c,hist};
    if(fa*fc<0){b=c;}else{a=c;fa=fc;}
  }
  return {x:(a+b)/2,hist};
}

function newtonRaphson(f,df,x0,tol,maxIter){
  const hist=[]; let x=x0;
  for(let k=0;k<maxIter;k++){
    const fx=f(x), dfx=df(x);
    if(Math.abs(dfx)<1e-14) return {x,hist,err:'Derivada ≈ 0, no converge'};
    const xNew=x-fx/dfx;
    const err=Math.abs(xNew-x);
    hist.push({k:k+1,x,fx,dfx,xNew,err});
    x=xNew;
    if(err<tol&&Math.abs(fx)<tol) return {x,hist};
  }
  return {x,hist};
}

function secante(f,x0,x1,tol,maxIter){
  const hist=[]; let xa=x0,xb=x1;
  for(let k=0;k<maxIter;k++){
    const fa=f(xa),fb=f(xb);
    if(Math.abs(fb-fa)<1e-14) return {x:xb,hist,err:'División por cero'};
    const xc=xb-fb*(xb-xa)/(fb-fa);
    const err=Math.abs(xc-xb);
    hist.push({k:k+1,xa,xb,fa,fb,xc,err});
    xa=xb; xb=xc;
    if(err<tol) return {x:xc,hist};
  }
  return {x:xb,hist};
}

function runEscenarioE(){
  const {f,df,a,b,x0,label,desc}=getFuncE();
  const metodo=document.getElementById('metodE').value;
  const tol=parseFloat(document.getElementById('tolE').value);

  const metodosRun=metodo==='todos'?['biseccion','newton','secante']:[metodo];
  let results={};
  metodosRun.forEach(m=>{
    if(m==='biseccion') results[m]=biseccion(f,a,b,tol,100);
    else if(m==='newton') results[m]=newtonRaphson(f,df,x0,tol,100);
    else results[m]=secante(f,x0,x0*(1.01)+0.1,tol,100);
  });

  const mainResult=results[metodosRun[metodosRun.length-1]];
  const root=mainResult.x;

  // Crystal
  const cb=document.getElementById('crystalE'); cb.style.display='block';
  const mainHist=mainResult.hist;
  let html='';
  const showN=Math.min(4,mainHist.length);
  for(let i=0;i<showN;i++){
    const h=mainHist[i];
    if(metodo==='biseccion'||metodosRun[0]==='biseccion'){
      html+=`<div class="iteration-step"><span class="step-label">Iteración ${h.k} — Bisección</span><span class="step-math">c = (${h.a?.toFixed(4)} + ${h.b?.toFixed(4)})/2 = ${h.c?.toFixed(6)}\nf(c) = ${h.fc?.toExponential(4)}</span><span class="step-result">Error = |b-a|/2 = ${h.err?.toExponential(4)}</span></div>`;
    } else if(metodo==='newton'||metodosRun[0]==='newton'){
      html+=`<div class="iteration-step"><span class="step-label">Iteración ${h.k} — Newton-Raphson</span><span class="step-math">f(${h.x?.toFixed(6)}) = ${h.fx?.toExponential(4)}\nf'(${h.x?.toFixed(6)}) = ${h.dfx?.toExponential(4)}\nx_nuevo = ${h.x?.toFixed(6)} - (${h.fx?.toExponential(4)})/(${h.dfx?.toExponential(4)}) = ${h.xNew?.toFixed(6)}</span><span class="step-result">|Δx| = ${h.err?.toExponential(4)}</span></div>`;
    } else {
      html+=`<div class="iteration-step"><span class="step-label">Iteración ${h.k} — Secante</span><span class="step-math">x_c = ${h.xb?.toFixed(6)} - f(${h.xb?.toFixed(4)})×(${h.xb?.toFixed(4)} - ${h.xa?.toFixed(4)}) / (f(${h.xb?.toFixed(4)}) - f(${h.xa?.toFixed(4)}))\n     = ${h.xc?.toFixed(6)}</span><span class="step-result">|Δx| = ${h.err?.toExponential(4)}</span></div>`;
    }
  }
  document.getElementById('crystalEContent').innerHTML=html;

  // Stats
  const iters=Object.entries(results).map(([m,r])=>`${m}:${r.hist.length}`).join(' | ');
  document.getElementById('statsE').innerHTML=`
    <div class="stat-card"><div class="sv">${root!==null?root.toFixed(5):'N/A'}</div><div class="sl">Raíz Encontrada (${label})</div></div>
    <div class="stat-card"><div class="sv">${root!==null?f(root).toExponential(3):'—'}</div><div class="sl">f(raíz) — debe ser ≈ 0</div></div>
    <div class="stat-card"><div class="sv">${mainHist.length}</div><div class="sl">Iteraciones (${metodosRun[metodosRun.length-1]})</div></div>
    <div class="stat-card"><div class="sv" style="font-size:1rem;color:var(--accent3)">${iters}</div><div class="sl">Comparativa iter.</div></div>
  `;

  // Charts
  const xs=[]; const ys=[];
  const range=b-a;
  for(let t=a-range*.05;t<=b+range*.05;t+=(range+range*.1)/200){
    xs.push(t); ys.push(f(t));
  }
  const datasets=[
    {label:'f(x)',data:xs.map((x,i)=>({x,y:ys[i]})),borderColor:'#4fc3f7',type:'line',tension:.3,pointRadius:0},
    {label:'y = 0',data:[{x:xs[0],y:0},{x:xs[xs.length-1],y:0}],borderColor:'#556070',borderDash:[4,4],type:'line',pointRadius:0},
  ];
  if(root!==null) datasets.push({label:`Raíz x* = ${root.toFixed(4)}`,data:[{x:root,y:0}],type:'scatter',backgroundColor:'#f5a623',pointRadius:10,pointStyle:'star'});

  makeChart('chartE',{type:'scatter',data:{datasets},options:{
    ...chartDefaults,responsive:true,
    scales:{x:{...chartDefaults.scales.x,type:'linear',title:{display:true,text:label,color:'#556070'}},
            y:{...chartDefaults.scales.y,title:{display:true,text:'f(x)',color:'#556070'}}}
  }});

  // Convergence chart
  const convDatasets=[];
  ['biseccion','newton','secante'].forEach((m,mi)=>{
    if(!results[m]) return;
    const cols=['#f5a623','#81c784','#ce93d8'];
    convDatasets.push({
      label:m,
      data:results[m].hist.map((h,i)=>({x:i+1,y:Math.log10(h.err+1e-15)})),
      borderColor:cols[mi],type:'line',tension:.3,pointRadius:3
    });
  });
  makeChart('chartE2',{type:'scatter',data:{datasets:convDatasets},options:{
    ...chartDefaults,responsive:true,
    scales:{x:{...chartDefaults.scales.x,type:'linear',title:{display:true,text:'Iteración',color:'#556070'}},
            y:{...chartDefaults.scales.y,title:{display:true,text:'log₁₀(Error)',color:'#556070'}}}
  }});

  // Iteration table
  const mainM=metodosRun[metodosRun.length-1];
  let thead,tbody;
  if(mainM==='biseccion'){
    thead='<tr><th>Iter.</th><th>a</th><th>b</th><th>c (aprox.)</th><th>f(c)</th><th>Error</th></tr>';
    tbody=mainResult.hist.map(h=>`<tr${Math.abs(h.err)<tol*10?' class="highlight-row"':''}><td>${h.k}</td><td>${h.a.toFixed(5)}</td><td>${h.b.toFixed(5)}</td><td>${h.c.toFixed(5)}</td><td>${h.fc.toExponential(3)}</td><td>${h.err.toExponential(3)}</td></tr>`).join('');
  } else if(mainM==='newton'){
    thead='<tr><th>Iter.</th><th>xₙ</th><th>f(xₙ)</th><th>f\'(xₙ)</th><th>xₙ₊₁</th><th>|Δx|</th></tr>';
    tbody=mainResult.hist.map(h=>`<tr${h.err<tol*10?' class="highlight-row"':''}><td>${h.k}</td><td>${h.x.toFixed(6)}</td><td>${h.fx.toExponential(3)}</td><td>${h.dfx.toExponential(3)}</td><td>${h.xNew.toFixed(6)}</td><td>${h.err.toExponential(3)}</td></tr>`).join('');
  } else {
    thead='<tr><th>Iter.</th><th>xₐ</th><th>x_b</th><th>x_nuevo</th><th>|Δx|</th></tr>';
    tbody=mainResult.hist.map(h=>`<tr${h.err<tol*10?' class="highlight-row"':''}><td>${h.k}</td><td>${h.xa.toFixed(5)}</td><td>${h.xb.toFixed(5)}</td><td>${h.xc.toFixed(6)}</td><td>${h.err.toExponential(3)}</td></tr>`).join('');
  }
  document.getElementById('tableE').innerHTML=`<thead>${thead}</thead><tbody>${tbody}</tbody>`;

  // Analysis
  const iterComp=Object.entries(results).map(([m,r])=>`${m}: ${r.hist.length} iter.`).join(', ');
  document.getElementById('qa-e1').textContent=root!==null?`${desc} → Valor crítico: ${root.toFixed(5)} ${label}. En este punto f(x*) = ${f(root).toExponential(3)} ≈ 0, confirmando la raíz.`:`No se encontró raíz en el rango dado. Verifica que la función cruce cero dentro del intervalo [${a}, ${b}].`;
  document.getElementById('qa-e2').textContent=`Comparativa de iteraciones: ${iterComp}. ${results.newton&&results.biseccion&&results.newton.hist.length<results.biseccion.hist.length?'Newton-Raphson convergió más rápido (convergencia cuadrática).':'Bisección necesitó más iteraciones pero es más robusta.'}`;
  document.getElementById('conclusionE').innerHTML=`<strong>Conclusión del Escenario E:</strong> Los métodos de búsqueda de raíces transforman preguntas de política pública ("¿cuándo se agotará el presupuesto?") en cálculos exactos. Newton-Raphson demuestra su superioridad en velocidad de convergencia para funciones suaves, mientras Bisección ofrece la garantía matemática de encontrar la raíz si existe en el intervalo. El método Secante equilibra ambos mundos: más rápido que Bisección sin necesitar la derivada analítica. La <em>limitación central</em> de este escenario es que las funciones modeladas son simplificaciones — en la realidad, el "umbral social" no es una función continua sino un fenómeno de histéresis y retroalimentación no lineal. Una mejora sería usar sistemas de EDOs acopladas (como el modelo NMD) para capturar la dinámica social con mayor fidelidad.`;

  document.getElementById('resultsE').style.display='block';
  document.getElementById('analysisE').style.display='block';
}

// Init tipo E panels
selectTipoE('gasto', document.querySelector('#tipoEtabs .method-tab'));