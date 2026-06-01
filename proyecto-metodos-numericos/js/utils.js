// ━━━━━━━━━━ NAVIGATION ━━━━━━━━━━
const sections=['home','scenA','scenB','scenC','scenD','scenE'];
function show(id,el){
  sections.forEach(s=>document.getElementById(s).classList.remove('active'));
  document.getElementById(id).classList.add('active');
  document.querySelectorAll('.nav-links a').forEach(a=>a.classList.remove('active'));
  if(el)el.classList.add('active');
  window.scrollTo({top:0,behavior:'smooth'});
}

// nav for scenario cards
document.querySelectorAll('.sc-card').forEach(c=>{
  const map={a:'scenA',b:'scenB',c:'scenC',d:'scenD',e:'scenE'};
  for(const [k,v] of Object.entries(map)){
    if(c.classList.contains(k)) c.onclick=()=>show(v,null);
  }
});

// chart registry
const charts={};
function mkChart(id,cfg){
  if(charts[id])charts[id].destroy();
  charts[id]=new Chart(document.getElementById(id).getContext('2d'),cfg);
}

function answerCard(q,a){
  return`<div class="answer-card"><div class="q">${q}</div><div class="a">${a}</div></div>`;
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  ESCENARIO A — Sistemas Lineales
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
function getMatrix(){
  const A=[[+document.getElementById('a00').value,+document.getElementById('a01').value,+document.getElementById('a02').value],
            [+document.getElementById('a10').value,+document.getElementById('a11').value,+document.getElementById('a12').value],
            [+document.getElementById('a20').value,+document.getElementById('a21').value,+document.getElementById('a22').value]];
  const b=[+document.getElementById('b0').value,+document.getElementById('b1').value,+document.getElementById('b2').value];
  return{A,b};
}

function luDecomp(A){
  const n=A.length;
  const L=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
  const U=A.map(r=>[...r]);
  const P=Array.from({length:n},(_,i)=>{const r=new Array(n).fill(0);r[i]=1;return r;});
  for(let k=0;k<n;k++){
    let maxVal=Math.abs(U[k][k]),maxRow=k;
    for(let i=k+1;i<n;i++) if(Math.abs(U[i][k])>maxVal){maxVal=Math.abs(U[i][k]);maxRow=i;}
    if(maxRow!==k){[U[k],U[maxRow]]=[U[maxRow],U[k]];[P[k],P[maxRow]]=[P[maxRow],P[k]];}
    for(let i=k+1;i<n;i++){
      L[i][k]=U[i][k]/U[k][k];
      for(let j=k;j<n;j++) U[i][j]-=L[i][k]*U[k][j];
    }
  }
  return{L,U,P};
}

function fwdSub(L,b){
  const n=L.length,y=new Array(n).fill(0);
  for(let i=0;i<n;i++){y[i]=b[i];for(let j=0;j<i;j++)y[i]-=L[i][j]*y[j];}
  return y;
}
function bkSub(U,y){
  const n=U.length,x=new Array(n).fill(0);
  for(let i=n-1;i>=0;i--){x[i]=y[i];for(let j=i+1;j<n;j++)x[i]-=U[i][j]*x[j];x[i]/=U[i][i];}
  return x;
}

function gaussianSolve(A,b){
  const n=A.length;
  const Aug=A.map((r,i)=>[...r,b[i]]);
  for(let k=0;k<n;k++){
    let maxR=k;
    for(let i=k+1;i<n;i++) if(Math.abs(Aug[i][k])>Math.abs(Aug[maxR][k]))maxR=i;
    [Aug[k],Aug[maxR]]=[Aug[maxR],Aug[k]];
    for(let i=k+1;i<n;i++){
      const f=Aug[i][k]/Aug[k][k];
      for(let j=k;j<=n;j++) Aug[i][j]-=f*Aug[k][j];
    }
  }
  const x=new Array(n).fill(0);
  for(let i=n-1;i>=0;i--){x[i]=Aug[i][n];for(let j=i+1;j<n;j++)x[i]-=Aug[i][j]*x[j];x[i]/=Aug[i][i];}
  return{x,iters:[{iter:1,x0:x[0].toFixed(4),x1:x[1].toFixed(4),x2:x[2].toFixed(4),err:'—'}]};
}

function jacobiSolve(A,b,tol,maxIter){
  const n=A.length;let x=new Array(n).fill(0);const iters=[];
  for(let it=1;it<=maxIter;it++){
    const xn=new Array(n).fill(0);
    for(let i=0;i<n;i++){let s=b[i];for(let j=0;j<n;j++) if(i!==j) s-=A[i][j]*x[j];xn[i]=s/A[i][i];}
    const err=Math.max(...xn.map((v,i)=>Math.abs(v-x[i])));
    iters.push({iter:it,x0:xn[0].toFixed(4),x1:xn[1].toFixed(4),x2:xn[2].toFixed(4),err:err.toExponential(3)});
    x=xn;if(err<tol)break;
  }
  return{x,iters};
}

function gaussSeidelSolve(A,b,tol,maxIter){
  const n=A.length;let x=new Array(n).fill(0);const iters=[];
  for(let it=1;it<=maxIter;it++){
    const xo=[...x];
    for(let i=0;i<n;i++){let s=b[i];for(let j=0;j<n;j++) if(i!==j) s-=A[i][j]*x[j];x[i]=s/A[i][i];}
    const err=Math.max(...x.map((v,i)=>Math.abs(v-xo[i])));
    iters.push({iter:it,x0:x[0].toFixed(4),x1:x[1].toFixed(4),x2:x[2].toFixed(4),err:err.toExponential(3)});
    if(err<tol)break;
  }
  return{x,iters};
}

function sorSolve(A,b,omega,tol,maxIter){
  const n=A.length;let x=new Array(n).fill(0);const iters=[];
  for(let it=1;it<=maxIter;it++){
    const xo=[...x];
    for(let i=0;i<n;i++){let s=b[i];for(let j=0;j<n;j++) if(i!==j) s-=A[i][j]*x[j];x[i]=(1-omega)*xo[i]+omega*s/A[i][i];}
    const err=Math.max(...x.map((v,i)=>Math.abs(v-xo[i])));
    iters.push({iter:it,x0:x[0].toFixed(4),x1:x[1].toFixed(4),x2:x[2].toFixed(4),err:err.toExponential(3)});
    if(err<tol)break;
  }
  return{x,iters};
}

function condNumber(A){
  const n=A.length;
  const normA=Math.sqrt(A.flat().reduce((s,v)=>s+v*v,0));
  const Aug=A.map((r,i)=>{const e=new Array(n).fill(0);e[i]=1;return[...r,...e];});
  for(let k=0;k<n;k++){
    let maxR=k;for(let i=k+1;i<n;i++) if(Math.abs(Aug[i][k])>Math.abs(Aug[maxR][k]))maxR=i;
    [Aug[k],Aug[maxR]]=[Aug[maxR],Aug[k]];
    const piv=Aug[k][k];for(let j=0;j<2*n;j++)Aug[k][j]/=piv;
    for(let i=0;i<n;i++){if(i===k)continue;const f=Aug[i][k];for(let j=0;j<2*n;j++)Aug[i][j]-=f*Aug[k][j];}
  }
  const Ainv=Aug.map(r=>r.slice(n));
  const normInv=Math.sqrt(Ainv.flat().reduce((s,v)=>s+v*v,0));
  return normA*normInv;
}

function solveA(){
  const{A,b}=getMatrix();
  const met=document.getElementById('metA').value;
  const tol=+document.getElementById('tolA').value;
  const maxI=+document.getElementById('maxiterA').value;
  const omega=+document.getElementById('omegaA').value;
  let res;
  if(met==='gauss') res=gaussianSolve(A,b);
  else if(met==='jacobi') res=jacobiSolve(A,b,tol,maxI);
  else if(met==='gauss_seidel') res=gaussSeidelSolve(A,b,tol,maxI);
  else res=sorSolve(A,b,omega,tol,maxI);

  const{x,iters}=res;
  const cond=condNumber(A);
  const zonas=['Norte','Centro','Sur'];
  document.getElementById('resultsA').style.display='block';

  document.getElementById('solA_vals').innerHTML=zonas.map((z,i)=>
    `<div class="result-box"><div class="rb-val">${x[i].toFixed(2)} L</div><div class="rb-lbl">Zona ${z}</div></div>`
  ).join('');

  const tbl=document.getElementById('tblA');
  tbl.innerHTML=`<thead><tr><th>Iter</th><th>x₁ (Norte)</th><th>x₂ (Centro)</th><th>x₃ (Sur)</th><th>Error</th></tr></thead><tbody>
    ${iters.map(r=>`<tr><td>${r.iter}</td><td>${r.x0}</td><td>${r.x1}</td><td>${r.x2}</td><td>${r.err}</td></tr>`).join('')}
  </tbody>`;

  const errs=iters.map(r=>r.err==='—'?null:+r.err.replace(/[^e\d\.\-]/g,''));
  mkChart('chartA',{
    type:'line',
    data:{labels:iters.map(r=>r.iter),datasets:[{label:'Error absoluto',data:errs.filter(v=>v!==null),borderColor:'#4361ee',backgroundColor:'rgba(67,97,238,.1)',fill:true,tension:.4,pointRadius:4}]},
    options:{responsive:true,maintainAspectRatio:false,scales:{y:{type:'logarithmic',title:{display:true,text:'Error (log)'}},x:{title:{display:true,text:'Iteración'}}},plugins:{legend:{position:'top'}}}
  });

  mkChart('chartA2',{
    type:'bar',
    data:{labels:['Zona Norte','Zona Centro','Zona Sur'],datasets:[
      {label:'Enviado (L)',data:x.map(v=>+v.toFixed(2)),backgroundColor:['rgba(67,97,238,.8)','rgba(114,9,183,.8)','rgba(247,37,133,.8)'],borderRadius:8}
    ]},
    options:{responsive:true,maintainAspectRatio:false,plugins:{legend:{position:'top'}},scales:{y:{title:{display:true,text:'Litros enviados'}}}}
  });

  const bloqueada=x.findIndex(v=>v<1);
  const minZona=zonas[x.indexOf(Math.min(...x))];
  const maxZona=zonas[x.indexOf(Math.max(...x))];
  const estable=cond<100;
  document.getElementById('answersA').innerHTML=[
    answerCard('¿Cuánto debe enviarse a cada zona?',`Norte: <strong>${x[0].toFixed(1)} L</strong>, Centro: <strong>${x[1].toFixed(1)} L</strong>, Sur: <strong>${x[2].toFixed(1)} L</strong>`),
    answerCard('¿Qué pasa si una ruta se bloquea?','Pon un coeficiente = 0 en la matriz y re-ejecuta. La demanda se redistribuye entre las rutas restantes, lo que puede sobrecargar otras.'),
    answerCard('¿Qué zona queda más afectada?',`La zona <strong>${minZona}</strong> recibe la menor cantidad (${Math.min(...x).toFixed(1)} L) y sería la más vulnerable ante un bloqueo.`),
    answerCard('¿El sistema es estable?',`Número de condición κ ≈ <strong>${cond.toFixed(1)}</strong>. ${estable?'Sistema <strong>bien condicionado</strong>: pequeños cambios producen pequeñas variaciones.':'Sistema <strong>mal condicionado</strong>: pequeños cambios en demanda generan grandes variaciones en la solución.'}`),
    answerCard('¿La solución cambia si la demanda sube?',`Sí. Con κ = ${cond.toFixed(1)}, un aumento del 5% en demanda puede provocar un error relativo de hasta ${(cond*0.05).toFixed(1)}% en la solución.`),
  ].join('');

  const metNombre={'gauss':'Eliminación Gaussiana + LU','jacobi':'Jacobi','gauss_seidel':'Gauss-Seidel','sor':'SOR'};
  document.getElementById('explainA').innerHTML=`
    <strong>¿Por qué este resultado?</strong><br>
    El método <strong>${metNombre[met]}</strong> resuelve el sistema Ax=b que modela las restricciones de distribución. 
    Cada ecuación representa una zona: la suma ponderada de lo enviado desde cada planta debe satisfacer la demanda exacta. 
    ${iters.length>1?`El método convergió en <strong>${iters.length} iteraciones</strong>, reduciendo el error hasta la tolerancia requerida.`:'Se obtuvo solución exacta mediante eliminación directa.'}
    El número de condición κ = ${cond.toFixed(1)} indica que el sistema es ${estable?'<strong>estable</strong> y confiable':'<strong>sensible</strong> a perturbaciones'}. 
    En contextos de crisis, un sistema bien condicionado significa que aunque la demanda varíe ligeramente, la distribución no colapsará.
  `;
  document.getElementById('resultsA').scrollIntoView({behavior:'smooth'});
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  ESCENARIO B — EDO
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
function solveB(){
  const R0=+document.getElementById('R0').value;
  const entrada=+document.getElementById('entrada').value;
  const consumo=+document.getElementById('consumoBase').value;
  const panico=+document.getElementById('panico').value;
  const dias=+document.getElementById('diasB').value;
  const h=+document.getElementById('hB').value;
  const crit=+document.getElementById('critico').value;
  const met=document.getElementById('metB').value;

  const f=(t,R)=>entrada - consumo*panico*(1+0.02*t);

  const N=Math.ceil(dias/h);
  const ts=[],euler=[],heun=[],rk4=[];
  let Re=R0,Rh=R0,Rr=R0;

  for(let i=0;i<=N;i++){
    const t=i*h;
    ts.push(+t.toFixed(3));
    euler.push(+Re.toFixed(2));
    heun.push(+Rh.toFixed(2));
    rk4.push(+Rr.toFixed(2));
    if(i<N){
      Re+=h*f(t,Re);
      const k1h=f(t,Rh),kp=Rh+h*k1h,k2h=f(t+h,kp);
      Rh+=h*(k1h+k2h)/2;
      const k1=f(t,Rr),k2=f(t+h/2,Rr+h*k1/2),k3=f(t+h/2,Rr+h*k2/2),k4=f(t+h,Rr+h*k3);
      Rr+=h*(k1+2*k2+2*k3+k4)/6;
    }
  }

  const critDayRK4=rk4.findIndex(v=>v<=crit);
  const critDayE=euler.findIndex(v=>v<=crit);

  document.getElementById('resultsB').style.display='block';

  const datasets=[];
  const show=met==='all';
  if(met==='euler'||show) datasets.push({label:'Euler',data:euler,borderColor:'#f72585',backgroundColor:'rgba(247,37,133,.07)',fill:false,tension:.3,pointRadius:0,borderWidth:2});
  if(met==='heun'||show) datasets.push({label:'Heun',data:heun,borderColor:'#4cc9f0',backgroundColor:'rgba(76,201,240,.07)',fill:false,tension:.3,pointRadius:0,borderWidth:2});
  if(met==='rk4'||show) datasets.push({label:'RK4',data:rk4,borderColor:'#4361ee',backgroundColor:'rgba(67,97,238,.1)',fill:true,tension:.3,pointRadius:0,borderWidth:2.5});

  const critLine={label:`Nivel crítico (${crit.toLocaleString()} L)`,data:ts.map(()=>crit),borderColor:'#e63946',borderDash:[6,4],pointRadius:0,borderWidth:2,fill:false};
  datasets.push(critLine);

  mkChart('chartB',{
    type:'line',
    data:{labels:ts,datasets},
    options:{responsive:true,maintainAspectRatio:false,scales:{y:{title:{display:true,text:'Reserva (litros)'},min:0},x:{title:{display:true,text:'Días'}}},plugins:{legend:{position:'top'}}}
  });

  const diffs=euler.map((v,i)=>Math.abs(v-rk4[i]));
  mkChart('chartB2',{
    type:'line',
    data:{labels:ts,datasets:[
      {label:'|Euler − RK4|',data:diffs,borderColor:'#f72585',backgroundColor:'rgba(247,37,133,.1)',fill:true,tension:.3,pointRadius:0},
      {label:'|Heun − RK4|',data:heun.map((v,i)=>Math.abs(v-rk4[i])),borderColor:'#4cc9f0',backgroundColor:'rgba(76,201,240,.1)',fill:true,tension:.3,pointRadius:0}
    ]},
    options:{responsive:true,maintainAspectRatio:false,scales:{y:{title:{display:true,text:'Diferencia absoluta'}},x:{title:{display:true,text:'Días'}}},plugins:{legend:{position:'top'}}}
  });

  const rows=rk4.filter((_,i)=>i%Math.max(1,Math.round(3/h))===0).map((v,i)=>{
    const d=(i*3).toFixed(1);
    return`<tr><td>${d}</td><td>${v.toLocaleString()}</td><td>${euler[Math.min(i*Math.round(3/h),euler.length-1)].toLocaleString()}</td><td>${(f(i*3,v)).toFixed(1)}</td><td>${v<=crit?'⚠️ CRÍTICO':'✅ OK'}</td></tr>`;
  });
  document.getElementById('tblB').innerHTML=`<thead><tr><th>Día</th><th>Reserva RK4 (L)</th><th>Reserva Euler (L)</th><th>dR/dt</th><th>Estado</th></tr></thead><tbody>${rows.join('')}</tbody>`;

  const diasCrit=critDayRK4>=0?(critDayRK4*h).toFixed(1):'No alcanza nivel crítico en el período';
  document.getElementById('answersB').innerHTML=[
    answerCard('¿En cuántos días la reserva llega a nivel crítico?',critDayRK4>=0?`En <strong>${diasCrit} días</strong> la reserva cae a ${crit.toLocaleString()} L según RK4.`:'La reserva <strong>no alcanza</strong> el nivel crítico en el período simulado.'),
    answerCard('¿Qué pasa si aumenta el consumo diario?','Aumenta el factor de pánico (≥ 1.5) y re-ejecuta. El día crítico se adelanta significativamente porque el término de consumo domina la EDO.'),
    answerCard('¿Qué pasa si se reduce el abastecimiento?','Disminuye "Entrada constante". Con entrada = 0 (bloqueo total) la reserva cae en caída libre siguiendo una curva casi lineal decreciente.'),
    answerCard('¿Qué método da aproximación más estable?','<strong>RK4</strong> es el más preciso: usa 4 evaluaciones por paso y es de orden 4. Heun (orden 2) mejora a Euler pero acumula más error para pasos grandes.'),
    answerCard('Diferencia Euler vs Heun vs RK4','Euler: 1 evaluación/paso, error O(h). Heun: 2 evaluaciones, error O(h²). RK4: 4 evaluaciones, error O(h⁴). Para h=1 día, la diferencia puede ser de cientos de litros.')
  ].join('');

  document.getElementById('explainB').innerHTML=`
    <strong>¿Por qué este resultado?</strong><br>
    La EDO <strong>R'(t) = ${entrada} − ${consumo}·${panico}·(1+0.02t)</strong> modela la reserva donde el consumo crece un 2%/día 
    por el efecto de pánico en la compra. Con entrada = ${entrada} L/día y consumo base = ${consumo} L/día × ${panico} (factor pánico), 
    el balance neto es <strong>${(entrada-consumo*panico).toFixed(0)} L/día</strong> 
    ${entrada<consumo*panico?'(negativo: la reserva <strong>disminuye</strong>)':'(positivo: la reserva <strong>aumenta</strong>)'}.<br><br>
    RK4 es el estándar de oro para EDOs por su excelente balance entre precisión y costo computacional. 
    La diferencia entre Euler y RK4 ilustra cómo un método de baja precisión puede subestimar o sobreestimar el día de crisis, 
    con consecuencias graves para la toma de decisiones.
  `;
  document.getElementById('resultsB').scrollIntoView({behavior:'smooth'});
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  ESCENARIO C — Interpolación
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
const defaultPtsC=[{x:1,y:8},{x:5,y:10},{x:10,y:13},{x:15,y:16},{x:20,y:19},{x:30,y:22}];
let ptsC=[...defaultPtsC.map(p=>({...p}))];

function renderPtsC(){
  document.getElementById('ptsC').innerHTML=`
    <div class="tbl-wrap"><table>
      <thead><tr><th>Día</th><th>Precio (Bs)</th><th></th></tr></thead>
      <tbody>${ptsC.map((p,i)=>`<tr>
        <td><input type="number" value="${p.x}" min="1" max="30" onchange="ptsC[${i}].x=+this.value" style="width:80px"></td>
        <td><input type="number" value="${p.y}" min="0" step="0.5" onchange="ptsC[${i}].y=+this.value" style="width:100px"></td>
        <td><button onclick="ptsC.splice(${i},1);renderPtsC()" style="background:none;border:none;cursor:pointer;color:#f72585;font-size:1.1rem">✕</button></td>
      </tr>`).join('')}
      </tbody>
    </table></div>`;
}
function addPtC(){ptsC.push({x:25,y:20.5});renderPtsC();}
function resetPtsC(){ptsC=defaultPtsC.map(p=>({...p}));renderPtsC();}
renderPtsC();

function lagrange(pts,xval){
  const n=pts.length;let y=0;
  for(let i=0;i<n;i++){
    let L=1;for(let j=0;j<n;j++) if(i!==j) L*=(xval-pts[j].x)/(pts[i].x-pts[j].x);
    y+=pts[i].y*L;
  }
  return y;
}

function newtonDivDiff(pts){
  const n=pts.length,x=pts.map(p=>p.x);
  const f=pts.map(p=>[p.y]);
  for(let j=1;j<n;j++) for(let i=0;i<n-j;i++) f[i].push((f[i+1][j-1]-f[i][j-1])/(x[i+j]-x[i]));
  return{x,coeffs:f[0].map((v,i)=>v)};
}
function newtonEval({x,coeffs},xval){
  const n=coeffs.length;let y=coeffs[n-1];
  for(let i=n-2;i>=0;i--) y=y*(xval-x[i])+coeffs[i];
  return y;
}

function splineCubic(pts){
  const n=pts.length,h=[];
  for(let i=0;i<n-1;i++) h.push(pts[i+1].x-pts[i].x);
  const A=Array.from({length:n},()=>new Array(n).fill(0));
  const rhs=new Array(n).fill(0);
  A[0][0]=1;A[n-1][n-1]=1;
  for(let i=1;i<n-1;i++){
    A[i][i-1]=h[i-1];A[i][i]=2*(h[i-1]+h[i]);A[i][i+1]=h[i];
    rhs[i]=3*((pts[i+1].y-pts[i].y)/h[i]-(pts[i].y-pts[i-1].y)/h[i-1]);
  }
  const c=[...rhs];const ab=A.map(r=>[...r]);
  for(let i=1;i<n;i++){const m=ab[i][i-1]/ab[i-1][i-1];ab[i][i]-=m*ab[i-1][i];c[i]-=m*c[i-1];}
  const sigma=new Array(n).fill(0);
  sigma[n-1]=c[n-1]/ab[n-1][n-1];
  for(let i=n-2;i>=0;i--) sigma[i]=(c[i]-ab[i][i+1]*sigma[i+1])/ab[i][i];
  return{pts,sigma,h};
}
function splineEval({pts,sigma,h},xval){
  const n=pts.length;
  let seg=n-2;
  for(let i=0;i<n-1;i++) if(xval>=pts[i].x&&xval<=pts[i+1].x){seg=i;break;}
  const dx=xval-pts[seg].x,hi=h[seg];
  const a=pts[seg].y,b=(pts[seg+1].y-pts[seg].y)/hi-hi*(2*sigma[seg]+sigma[seg+1])/3;
  const c=sigma[seg],d=(sigma[seg+1]-sigma[seg])/(3*hi);
  return a+b*dx+c*dx*dx+d*dx*dx*dx;
}

function solveC(){
  const sorted=[...ptsC].sort((a,b)=>a.x-b.x);
  const evalDay=+document.getElementById('evalDayC').value;
  const prod=document.getElementById('productoC').value;
  const newtonModel=newtonDivDiff(sorted);
  const splineModel=splineCubic(sorted);

  const days=[];for(let d=sorted[0].x;d<=sorted[sorted.length-1].x;d+=0.5)days.push(d);
  const lagVals=days.map(d=>lagrange(sorted,d));
  const newVals=days.map(d=>newtonEval(newtonModel,d));
  const splVals=days.map(d=>splineEval(splineModel,d));

  const lagEval=lagrange(sorted,evalDay);
  const newEval=newtonEval(newtonModel,evalDay);
  const splEval=splineEval(splineModel,evalDay);

  document.getElementById('resultsC').style.display='block';

  mkChart('chartC',{
    type:'line',
    data:{
      labels:days,
      datasets:[
        {label:'Datos reales',data:sorted.map(p=>({x:p.x,y:p.y})),type:'scatter',borderColor:'#1a1f3c',backgroundColor:'#1a1f3c',pointRadius:8,pointHoverRadius:10,showLine:false},
        {label:'Lagrange',data:lagVals,borderColor:'#4361ee',backgroundColor:'transparent',tension:.3,pointRadius:0,borderWidth:2},
        {label:'Newton',data:newVals,borderColor:'#f72585',backgroundColor:'transparent',tension:.3,pointRadius:0,borderWidth:2,borderDash:[5,3]},
        {label:'Splines cúbicos',data:splVals,borderColor:'#06d6a0',backgroundColor:'rgba(6,214,160,.08)',tension:.4,pointRadius:0,borderWidth:2.5,fill:true},
        {label:`Día ${evalDay} evaluado`,data:days.map(d=>Math.abs(d-evalDay)<0.26?splEval:null),type:'scatter',borderColor:'#ffd60a',backgroundColor:'#ffd60a',pointRadius:10,pointStyle:'triangle',showLine:false},
      ]
    },
    options:{responsive:true,maintainAspectRatio:false,scales:{y:{title:{display:true,text:`Precio ${prod} (Bs)`}},x:{title:{display:true,text:'Día del mes'}}},plugins:{legend:{position:'top'}}}
  });

  mkChart('chartC2',{
    type:'bar',
    data:{labels:['Lagrange','Newton','Splines'],datasets:[
      {label:`Precio en día ${evalDay} (Bs)`,data:[lagEval,newEval,splEval].map(v=>+v.toFixed(3)),backgroundColor:['rgba(67,97,238,.8)','rgba(247,37,133,.8)','rgba(6,214,160,.8)'],borderRadius:8}
    ]},
    options:{responsive:true,maintainAspectRatio:false,plugins:{legend:{position:'top'}},scales:{y:{title:{display:true,text:'Precio (Bs)'}}}}
  });

  const selDays=[1,5,7,10,12,15,18,20,25,30].filter(d=>d>=sorted[0].x&&d<=sorted[sorted.length-1].x);
  document.getElementById('tblC').innerHTML=`<thead><tr><th>Día</th><th>Lagrange (Bs)</th><th>Newton (Bs)</th><th>Splines (Bs)</th><th>Dato real</th></tr></thead>
    <tbody>${selDays.map(d=>{
      const real=sorted.find(p=>p.x===d);
      return`<tr><td>${d}</td><td>${lagrange(sorted,d).toFixed(3)}</td><td>${newtonEval(newtonModel,d).toFixed(3)}</td><td>${splineEval(splineModel,d).toFixed(3)}</td><td>${real?real.y+'✔':'—'}</td></tr>`;
    }).join('')}</tbody>`;

  const maxInc=sorted[sorted.length-1].y-sorted[0].y;
  const pct=((maxInc/sorted[0].y)*100).toFixed(1);
  const confiable=sorted.length>=5&&(sorted[sorted.length-1].x-sorted[0].x)<=30;

  document.getElementById('answersC').innerHTML=[
    answerCard('¿Cuál sería el precio en el día evaluado?',`Lagrange: <strong>${lagEval.toFixed(2)} Bs</strong> | Newton: <strong>${newEval.toFixed(2)} Bs</strong> | Splines: <strong>${splEval.toFixed(2)} Bs</strong>`),
    answerCard('¿Cómo se comporta la curva durante el mes?',`Tendencia <strong>alcista</strong>: de ${sorted[0].y} Bs a ${sorted[sorted.length-1].y} Bs, un incremento del <strong>${pct}%</strong> en el mes.`),
    answerCard('¿Cuál método es más confiable?',`Los <strong>Splines cúbicos</strong> son los más confiables: garantizan suavidad, evitan oscilaciones de Runge y son locales (cada segmento depende solo de sus vecinos).`),
    answerCard('¿Qué pasa si los datos son muy dispersos?',`Lagrange y Newton pueden oscilar violentamente entre puntos (fenómeno de Runge). Los Splines mantienen estabilidad local incluso con datos irregulares.`),
    answerCard('¿Qué producto tuvo mayor incremento?',`<strong>${prod}</strong>: subió <strong>${maxInc.toFixed(1)} Bs</strong> (${pct}%). Ajusta los datos para comparar con otros productos.`),
  ].join('');

  document.getElementById('explainC').innerHTML=`
    <strong>¿Por qué este resultado?</strong><br>
    Los tres métodos usan los mismos ${sorted.length} datos pero con estrategias distintas. 
    <strong>Lagrange</strong> construye un polinomio global de grado ${sorted.length-1}: preciso en los nodos pero susceptible a oscilaciones extremas en los bordes (fenómeno de Runge). 
    <strong>Newton con diferencias divididas</strong> es algebraicamente equivalente a Lagrange pero más eficiente para agregar nuevos puntos. 
    <strong>Splines cúbicos</strong> ajustan polinomios de grado 3 por tramos con condiciones de continuidad en derivadas primera y segunda, 
    lo que produce curvas suaves y estables, el método preferido en economía y análisis de precios. 
    La diferencia entre métodos en el día evaluado (≈${Math.abs(lagEval-splEval).toFixed(2)} Bs) representa la <em>incertidumbre de interpolación</em>.
  `;
  document.getElementById('resultsC').scrollIntoView({behavior:'smooth'});
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  ESCENARIO D — Integración
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
function buildPriceFunc(p0,pf,dias,type){
  if(type==='lineal') return t=>p0+(pf-p0)*t/dias;
  if(type==='exponencial') return t=>p0*Math.exp(Math.log(pf/p0)*t/dias);
  return t=>{
    if(t<10) return p0;
    if(t<20) return p0+(pf-p0)*0.4;
    return pf;
  };
}

function trapecio(f,a,b,n){
  const h=(b-a)/n;let s=(f(a)+f(b))/2;
  for(let i=1;i<n;i++) s+=f(a+i*h);
  return s*h;
}
function simpson13(f,a,b,n){
  if(n%2!==0)n++;
  const h=(b-a)/n;let s=f(a)+f(b);
  for(let i=1;i<n;i++) s+=(i%2===0?2:4)*f(a+i*h);
  return s*h/3;
}
function simpson38(f,a,b,n){
  if(n%3!==0)n+=3-n%3;
  const h=(b-a)/n;let s=f(a)+f(b);
  for(let i=1;i<n;i++) s+=(i%3===0?2:3)*f(a+i*h);
  return s*3*h/8;
}

function solveD(){
  const ingreso=+document.getElementById('ingresoD').value;
  const cant=+document.getElementById('cantD').value;
  const p0=+document.getElementById('p0D').value;
  const pf=+document.getElementById('pfD').value;
  const dias=+document.getElementById('diasD').value;
  const type=document.getElementById('growthD').value;
  const priceF=buildPriceFunc(p0,pf,dias,type);
  const gastoF=t=>priceF(t)*cant;
  const n=dias%2===0?dias:(dias-1);

  const gastTrap=trapecio(gastoF,0,dias,dias);
  const gastS13=simpson13(gastoF,0,dias,n);
  const gastS38=simpson38(gastoF,0,dias,n%3===0?n:n+(3-n%3));
  const gastConst=p0*cant*dias;

  const days=Array.from({length:dias+1},(_,i)=>i);
  const prices=days.map(d=>priceF(d));
  const cumGasto=days.map((_,i)=>i===0?0:trapecio(gastoF,0,i,Math.max(i,1)));

  document.getElementById('resultsD').style.display='block';

  mkChart('chartD',{
    type:'line',
    data:{
      labels:days,
      datasets:[
        {label:'Precio diario (Bs)',data:prices,borderColor:'#f72585',backgroundColor:'rgba(247,37,133,.08)',fill:true,tension:.3,yAxisID:'y1',pointRadius:3},
        {label:'Gasto acumulado (Bs)',data:cumGasto,borderColor:'#4361ee',backgroundColor:'rgba(67,97,238,.08)',fill:true,tension:.3,yAxisID:'y2',pointRadius:2},
        {label:'Ingreso mensual (Bs)',data:days.map(()=>ingreso),borderColor:'#06d6a0',borderDash:[6,4],pointRadius:0,yAxisID:'y2',borderWidth:2},
      ]
    },
    options:{responsive:true,maintainAspectRatio:false,
      scales:{
        y1:{type:'linear',position:'left',title:{display:true,text:'Precio Bs/día'}},
        y2:{type:'linear',position:'right',title:{display:true,text:'Gasto acumulado Bs'},grid:{drawOnChartArea:false}}
      },plugins:{legend:{position:'top'}}}
  });

  mkChart('chartD2',{
    type:'bar',
    data:{labels:['Trapecios','Simpson 1/3','Simpson 3/8','Sin inflación'],datasets:[
      {label:'Gasto mensual (Bs)',data:[gastTrap,gastS13,gastS38,gastConst].map(v=>+v.toFixed(2)),
       backgroundColor:['rgba(67,97,238,.8)','rgba(114,9,183,.8)','rgba(247,37,133,.8)','rgba(6,214,160,.8)'],borderRadius:8}
    ]},
    options:{responsive:true,maintainAspectRatio:false,plugins:{legend:{position:'top'}},scales:{y:{title:{display:true,text:'Bs'}}}}
  });

  const step=Math.max(1,Math.floor(dias/10));
  const rows=days.filter(d=>d%step===0||d===dias).map(d=>{
    const g=d===0?0:trapecio(gastoF,0,d,Math.max(d,1));
    return`<tr><td>${d}</td><td>${priceF(d).toFixed(2)}</td><td>${gastoF(d).toFixed(2)}</td><td>${g.toFixed(2)}</td><td>${g>ingreso?'⚠️ Supera':'✅ OK'}</td></tr>`;
  });
  document.getElementById('tblD').innerHTML=`<thead><tr><th>Día</th><th>Precio (Bs)</th><th>Gasto diario (Bs)</th><th>Acumulado (Bs)</th><th>vs Ingreso</th></tr></thead><tbody>${rows.join('')}</tbody>`;

  const perdida=gastTrap-gastConst;
  const pctPerd=(perdida/gastConst*100).toFixed(1);
  const diaSupera=cumGasto.findIndex(v=>v>ingreso);
  document.getElementById('summaryD').innerHTML=[
    `<div class="result-box"><div class="rb-val">${gastTrap.toFixed(0)} Bs</div><div class="rb-lbl">Gasto real del mes (Trapecios)</div></div>`,
    `<div class="result-box"><div class="rb-val">${gastConst.toFixed(0)} Bs</div><div class="rb-lbl">Gasto sin inflación</div></div>`,
    `<div class="result-box"><div class="rb-val">${perdida.toFixed(0)} Bs</div><div class="rb-lbl">Pérdida poder adquisitivo (${pctPerd}%)</div></div>`,
    diaSupera>0?`<div class="result-box" style="border-color:#f72585"><div class="rb-val" style="color:#f72585">Día ${diaSupera}</div><div class="rb-lbl">⚠️ Gasto supera ingreso</div></div>`:'',
  ].join('');

  const bestMethod=Math.abs(gastS13-gastS38)<Math.abs(gastTrap-gastS38)?'Simpson 1/3':'Trapecios';
  document.getElementById('answersD').innerHTML=[
    answerCard('¿Cuánto gastó la familia durante el mes?',`<strong>${gastTrap.toFixed(2)} Bs</strong> (Trapecios) | ${gastS13.toFixed(2)} Bs (Simpson 1/3) | ${gastS38.toFixed(2)} Bs (Simpson 3/8)`),
    answerCard('¿Cuánto hubiera gastado sin inflación?',`Con precio fijo de ${p0} Bs: <strong>${gastConst.toFixed(2)} Bs</strong>. La inflación costó <strong>${perdida.toFixed(2)} Bs</strong> adicionales.`),
    answerCard('¿Cuál fue la pérdida del poder adquisitivo?',`<strong>${pctPerd}%</strong> del gasto base. El presupuesto familiar se erosionó ${perdida.toFixed(0)} Bs por la suba de precios.`),
    answerCard('¿Qué método de integración fue más preciso?',`<strong>Simpson 1/3</strong> es generalmente más preciso que Trapecios (error O(h⁴) vs O(h²)). Para funciones suaves, la diferencia es mínima.`),
    answerCard('¿Qué producto afecta más al gasto?',`El producto analizado representa <strong>${((gastTrap/ingreso)*100).toFixed(1)}%</strong> del ingreso mensual. Cambia la cantidad diaria para comparar productos.`),
  ].join('');

  document.getElementById('explainD').innerHTML=`
    <strong>¿Por qué este resultado?</strong><br>
    El gasto acumulado es la integral ∫₀³⁰ P(t)·${cant} dt, donde P(t) es la curva de precios de tipo <strong>${type}</strong>. 
    <strong>Trapecios</strong> aproxima el área con trapecios de ancho h, con error O(h²). 
    <strong>Simpson 1/3</strong> usa parábolas por cada par de subintervalos, con error O(h⁴), generalmente más exacto. 
    <strong>Simpson 3/8</strong> usa cúbicas cada 3 subintervalos, eficiente cuando n es múltiplo de 3. 
    La diferencia entre métodos (${Math.abs(gastTrap-gastS13).toFixed(2)} Bs) es pequeña porque la curva de precios es suave; 
    con datos irregulares (escalones), Trapecios puede ser más robusto. 
    La pérdida de poder adquisitivo de <strong>${perdida.toFixed(0)} Bs</strong> representa el costo real de la crisis para una familia típica.
  `;
  document.getElementById('resultsD').scrollIntoView({behavior:'smooth'});
}

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
//  ESCENARIO E — Raíces
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
const defaultParamsE={
  gasto:{ingreso:2500,p0:8,pf:22,cant:2,dias:30,label:'Día t donde gasto(t) = ingreso',xmin:1,xmax:50},
  combustible:{R0:50000,entrada:3000,consumo:5000,crit:5000,label:'Tasa entrada k donde R(30)=crit',xmin:0,xmax:10000},
  social:{a:0.002,b:0.05,N0:900,M0:100,crit:0.5,label:'Parámetro a donde M(t)=N(t)',xmin:0.0001,xmax:0.01},
};

function updateEparams(){
  const tipo=document.getElementById('problemaE').value;
  const p=defaultParamsE[tipo];
  const fields={
    gasto:[{id:'eIngreso',label:'Ingreso familiar (Bs)',val:p.ingreso},{id:'eP0',label:'Precio inicial (Bs)',val:p.p0},{id:'ePf',label:'Precio final (Bs)',val:p.pf},{id:'eCant',label:'Cantidad diaria',val:p.cant}],
    combustible:[{id:'eR0',label:'Reserva inicial (L)',val:p.R0},{id:'eCons',label:'Consumo (L/día)',val:p.consumo},{id:'eCrit',label:'Nivel crítico (L)',val:p.crit},{id:'eDias',label:'Días',val:p.dias||30}],
    social:[{id:'eA',label:'Tasa contagio a',val:p.a},{id:'eB',label:'Recuperación b',val:p.b},{id:'eN0',label:'Neutrales N₀',val:p.N0},{id:'eM0',label:'Manifestantes M₀',val:p.M0}],
  };
  document.getElementById('paramsE').innerHTML=(fields[tipo]||[]).map(f=>
    `<div><label>${f.label}</label><input type="number" id="${f.id}" value="${f.val}" step="any"></div>`
  ).join('');
}
updateEparams();

function getFuncE(){
  const tipo=document.getElementById('problemaE').value;
  if(tipo==='gasto'){
    const ingreso=+(document.getElementById('eIngreso')?.value||2500);
    const p0=+(document.getElementById('eP0')?.value||8);
    const pf=+(document.getElementById('ePf')?.value||22);
    const cant=+(document.getElementById('eCant')?.value||2);
    const dias=30;
    const prF=t=>p0+(pf-p0)*t/dias;
    const f=t=>trapecio(s=>prF(s)*cant,0,t,Math.max(1,Math.round(t)))-ingreso;
    return{f,label:'Gasto acumulado − Ingreso',xmin:1,xmax:50,xunit:'días'};
  }
  if(tipo==='combustible'){
    const R0=+(document.getElementById('eR0')?.value||50000);
    const cons=+(document.getElementById('eCons')?.value||5000);
    const crit=+(document.getElementById('eCrit')?.value||5000);
    const f=k=>(R0+30*k-30*cons)-crit;
    return{f,label:'R(30,k) − Nivel crítico',xmin:0,xmax:10000,xunit:'L/día (tasa entrada)'};
  }
  const a=+(document.getElementById('eA')?.value||0.002);
  const b=+(document.getElementById('eB')?.value||0.05);
  const N0=+(document.getElementById('eN0')?.value||900);
  const M0=+(document.getElementById('eM0')?.value||100);
  function simSocial(aVal){
    let N=N0,M=M0,D=10;
    const h=0.5,c=0.03,k=0.05,r=0.1;
    for(let t=0;t<30;t+=h){
      const dN=-aVal*N*M+b*D, dM=aVal*N*M-c*M*D, dD=k*M-r*D;
      N+=h*dN;M+=h*dM;D+=h*dD;
    }
    return M-N;
  }
  return{f:simSocial,label:'M(30) − N(30) (manifestantes vs neutrales)',xmin:0.0001,xmax:0.005,xunit:'tasa de contagio a'};
}

function biseccion(f,a,b,tol,maxI){
  const iters=[];
  if(f(a)*f(b)>0) return{root:null,iters:[{iter:1,a,b,c:(a+b)/2,fc:'sin cruce',err:'—'}]};
  let lo=a,hi=b;
  for(let i=1;i<=maxI;i++){
    const c=(lo+hi)/2,fc=f(c),err=(hi-lo)/2;
    iters.push({iter:i,a:lo.toFixed(6),b:hi.toFixed(6),c:c.toFixed(6),fc:fc.toFixed(6),err:err.toExponential(3)});
    if(err<tol) return{root:c,iters};
    if(f(lo)*fc<0)hi=c;else lo=c;
  }
  return{root:(+iters[iters.length-1].c),iters};
}

function newtonR(f,x0,tol,maxI){
  const h=1e-7;const iters=[];let x=x0;
  for(let i=1;i<=maxI;i++){
    const fx=f(x),dfx=(f(x+h)-f(x-h))/(2*h);
    if(Math.abs(dfx)<1e-14) break;
    const xn=x-fx/dfx,err=Math.abs(xn-x);
    iters.push({iter:i,x:x.toFixed(6),fx:fx.toFixed(6),dfx:dfx.toFixed(6),xn:xn.toFixed(6),err:err.toExponential(3)});
    x=xn;if(err<tol)break;
  }
  return{root:x,iters};
}

function secante(f,x0,x1,tol,maxI){
  const iters=[];
  for(let i=1;i<=maxI;i++){
    const fx0=f(x0),fx1=f(x1),den=fx1-fx0;
    if(Math.abs(den)<1e-14) break;
    const xn=x1-fx1*(x1-x0)/den,err=Math.abs(xn-x1);
    iters.push({iter:i,x0:x0.toFixed(6),x1:x1.toFixed(6),xn:xn.toFixed(6),fxn:f(xn).toFixed(6),err:err.toExponential(3)});
    x0=x1;x1=xn;if(err<tol)break;
  }
  return{root:x1,iters};
}

function solveE(){
  const{f,label,xmin,xmax,xunit}=getFuncE();
  const tol=+document.getElementById('tolE').value;
  const maxI=+document.getElementById('maxitE').value;
  const met=document.getElementById('metE').value;

  const xmid=(xmin+xmax)/2;
  const resB=biseccion(f,xmin,xmax,tol,maxI);
  const resN=newtonR(f,xmid,tol,maxI);
  const resS=secante(f,xmin,xmax,tol,maxI);

  const root=resB.root??resN.root??resS.root??null;

  document.getElementById('resultsE').style.display='block';

  const xs=[],ys=[];
  const steps=200;
  for(let i=0;i<=steps;i++){const x=xmin+(xmax-xmin)*i/steps;xs.push(x);ys.push(f(x));}

  const datasets=[
    {label:`f(x) = ${label}`,data:ys,borderColor:'#4361ee',backgroundColor:'rgba(67,97,238,.08)',fill:true,tension:.3,pointRadius:0,borderWidth:2},
    {label:'f(x) = 0',data:xs.map(()=>0),borderColor:'#aaa',borderDash:[4,4],pointRadius:0,borderWidth:1},
  ];
  if(root!==null) datasets.push({label:`Raíz ≈ ${root.toFixed(4)} ${xunit}`,data:xs.map(x=>Math.abs(x-root)<(xmax-xmin)/steps*2?f(root):null),type:'scatter',backgroundColor:'#f72585',pointRadius:12,pointStyle:'star',showLine:false});

  mkChart('chartE',{
    type:'line',
    data:{labels:xs.map(v=>v.toFixed(3)),datasets},
    options:{responsive:true,maintainAspectRatio:false,scales:{y:{title:{display:true,text:'f(x)'}},x:{title:{display:true,text:xunit},ticks:{maxTicksLimit:12}}},plugins:{legend:{position:'top'}}}
  });

  mkChart('chartE2',{
    type:'line',
    data:{labels:Array.from({length:Math.max(resB.iters.length,resN.iters.length,resS.iters.length)},(_,i)=>i+1),
      datasets:[
        {label:`Bisección (${resB.iters.length} iter)`,data:resB.iters.map(r=>r.err==='—'?null:+r.err),borderColor:'#4361ee',fill:false,tension:.2,pointRadius:4},
        {label:`Newton-Raphson (${resN.iters.length} iter)`,data:resN.iters.map(r=>+r.err),borderColor:'#f72585',fill:false,tension:.2,pointRadius:4},
        {label:`Secante (${resS.iters.length} iter)`,data:resS.iters.map(r=>+r.err),borderColor:'#06d6a0',fill:false,tension:.2,pointRadius:4},
      ]
    },
    options:{responsive:true,maintainAspectRatio:false,scales:{y:{type:'logarithmic',title:{display:true,text:'Error (log)'}},x:{title:{display:true,text:'Iteración'}}},plugins:{legend:{position:'top'}}}
  });

  const bestRes=resN.iters.length<=resB.iters.length?resN:resB;
  document.getElementById('tblE').innerHTML=`
    <thead><tr><th>Iter</th><th>a / x₀</th><th>b / x₁</th><th>c / xₙ</th><th>f(c)</th><th>Error</th></tr></thead>
    <tbody>${[...resB.iters,...resN.iters.map(r=>({...r,a:r.x,b:r.dfx,c:r.xn,fc:r.fxn})),...resS.iters.map(r=>({...r,a:r.x0,b:r.x1,c:r.xn,fc:r.fxn}))].slice(0,20).map(r=>
    `<tr><td>${r.iter}</td><td>${r.a}</td><td>${r.b}</td><td>${r.c}</td><td>${r.fc}</td><td>${r.err}</td></tr>`).join('')}</tbody>`;

  const rootStr=root!==null?root.toFixed(4):'No encontrada';
  const tipo=document.getElementById('problemaE').value;
  const interpretation={
    gasto:`El gasto acumulado iguala al ingreso en el <strong>día ${rootStr}</strong>. A partir de ese momento la familia entra en déficit.`,
    combustible:`La tasa mínima de reposición para no llegar a nivel crítico es <strong>${rootStr} L/día</strong>. Por debajo de este valor, la planta colapsa.`,
    social:`Con tasa de contagio a = <strong>${rootStr}</strong> los manifestantes igualan a los neutrales: umbral de masificación del conflicto.`,
  };

  document.getElementById('answersE').innerHTML=[
    answerCard('¿Dónde está el umbral crítico?',interpretation[tipo]||`Raíz encontrada en x = <strong>${rootStr}</strong>`),
    answerCard('¿Qué método converge más rápido?',`Newton-Raphson: <strong>${resN.iters.length} iteraciones</strong>. Bisección: <strong>${resB.iters.length}</strong>. Secante: <strong>${resS.iters.length}</strong>. Newton es superlineal (orden ≈2), pero requiere derivada.`),
    answerCard('¿Cuál es más robusto?',`<strong>Bisección</strong> garantiza convergencia si hay un cruce de signo, aunque es más lento. Newton puede fallar si f\'(x)≈0 o la condición inicial está lejos de la raíz.`),
    answerCard('¿Qué orden de convergencia tiene cada uno?','Bisección: orden 1 (lineal). Newton-Raphson: orden 2 (cuadrático). Secante: orden ≈1.618 (superlineal). Mayor orden = menos iteraciones para la misma tolerancia.'),
    answerCard('¿Qué pasa si la condición inicial cambia?','Newton y Secante son sensibles al punto de inicio. Si f\'(x₀)≈0 o el punto está en una zona sin cruce, pueden divergir. Bisección solo necesita que f(a)·f(b)<0.'),
  ].join('');

  document.getElementById('explainE').innerHTML=`
    <strong>¿Por qué este resultado?</strong><br>
    Se busca la raíz de <strong>f(x) = ${label}</strong>, es decir, el valor de x donde el sistema "cambia de estado". 
    <strong>Bisección</strong> divide iterativamente el intervalo [${xmin}, ${xmax}] a la mitad, garantizando convergencia lineal. 
    Convergió en <strong>${resB.iters.length} iteraciones</strong>. 
    <strong>Newton-Raphson</strong> usa la tangente a la curva: xₙ₊₁ = xₙ − f(xₙ)/f\'(xₙ), con convergencia cuadrática. 
    Convergió en <strong>${resN.iters.length} iteraciones</strong>. 
    <strong>Secante</strong> aproxima la derivada con dos puntos previos, convergencia de orden φ≈1.618. 
    La raíz encontrada x ≈ <strong>${rootStr} ${xunit}</strong> es el umbral donde el sistema cambia de fase: 
    cruzar este punto significa que la crisis pasa de manejable a crítica.
  `;
  document.getElementById('resultsE').scrollIntoView({behavior:'smooth'});
}
