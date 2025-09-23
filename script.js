// --- CONFIGURAÇÃO & SETUP ---
const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d', { alpha: false });
// --- INÍCIO DA IMPLEMENTAÇÃO: Tooltip no Canvas ---
const canvasTooltip = document.getElementById('canvasTooltip');
// --- FIM DA IMPLEMENTAÇÃO ---

// --- PARÂMETROS FÍSICOS DO DOMÍNIO (SI) ---
const DOMAIN_WIDTH = 2.0; // Largura do domínio em metros
let DOMAIN_HEIGHT;       // Altura será calculada para manter o aspect ratio

let canvasW = Math.floor(window.innerWidth > 1200 ? Math.min(window.innerWidth - 420, 1280) : window.innerWidth - 40);
let canvasH = Math.floor(window.innerHeight - 80);
canvas.width = canvasW;
canvas.height = canvasH;

let NX = 180; // Resolução da grade em X
let NY;       // Resolução da grade em Y

// Células da grade (cálculo inicial)
let dx, dy; // Tamanho de cada célula da grade em metros

function updateGridDimensions() {
    DOMAIN_HEIGHT = DOMAIN_WIDTH * canvasH / canvasW;
    NY = Math.round(NX * DOMAIN_HEIGHT / DOMAIN_WIDTH);
    dx = DOMAIN_WIDTH / NX;
    dy = DOMAIN_HEIGHT / NY;
    document.getElementById('resLabel').textContent = `${NX} x ${NY}`;
}
updateGridDimensions();


// --- CAMPOS FÍSICOS ---
let size, u, v, u_prev, v_prev, r, g, b, r_prev, g_prev, b_prev, obstacle, div, p, vort;
let temp, temp_prev, density, visc_field, thermal_diff_field;

function alloc() {
  size = (NX + 2) * (NY + 2);
  u = new Float32Array(size); v = new Float32Array(size);
  u_prev = new Float32Array(size); v_prev = new Float32Array(size);
  r = new Float32Array(size); g = new Float32Array(size); b = new Float32Array(size);
  r_prev = new Float32Array(size); g_prev = new Float32Array(size); b_prev = new Float32Array(size);
  obstacle = new Uint8Array(size);
  div = new Float32Array(size); p = new Float32Array(size); vort = new Float32Array(size);
  
  temp = new Float32Array(size); temp_prev = new Float32Array(size);
  density = new Float32Array(size); // Densidade local ρ(T) para empuxo
  visc_field = new Float32Array(size); // Viscosidade dinâmica (μ) em Pa·s
  thermal_diff_field = new Float32Array(size); // Difusividade térmica (α) em m²/s
}

function IX(i, j) { return i + j * (NX + 2); }
alloc();

// --- PARÂMETROS DE SIMULAÇÃO & ESTADO ---
let paused = false;
let injectionVelocity = 5.0;  // Velocidade de injeção em m/s
let injectionTemp = 800;      // Temperatura de injeção em Kelvin (K)
let vorticityStrength = 0.15;
let injectionOn = true;
let injectionColor = { r: 255, g: 80, b: 0 };
let obstacles = [];
let draggedObstacleIndex = -1;
let vorticityVisScale = 10.0;
let pressureSensitivity = 1.0;
let temperatureSensitivity = 1.0;
setupSlider('tempSensRange', 'tempSensVal', v => temperatureSensitivity = v);
let showPressureContours = false;

// --- CORREÇÃO: Variáveis de gravação movidas para o topo para evitar ReferenceError ---
let isRecording = false;
let recordingDuration = 5.0;
let recordingFrequency = 20;
let simulationTimeRecorded = 0;
let recordedData = [];
let originalRecordedData = [];
let simulationParametersSnapshot = {};

// --- NOVO: Variável para armazenar o histórico de gravações ---
let recordingHistory = [];

// --- Flag para controlar a renderização dos gráficos de correlação ---
let correlationChartsRendered = false;


// --- CONSTANTES E PROPRIEDADES FÍSICAS (SI) ---
// Propriedades do fluido (padrão: ar)
let thermalConductivity_k = 0.0257; // Condutividade térmica [W/(m·K)]
let specificHeat_cp = 1005;         // Calor específico a pressão constante [J/(kg·K)]
let base_dynamic_viscosity = 1.81e-5; // Viscosidade dinâmica (μ) base para ar a 293K [Pa·s ou N·s/m²]
const R_specific_air = 287.058;     // Constante de gás específica para o ar [J/(kg·K)]

// Propriedades para Lei de Sutherland (viscosidade de gases)
const sutherland_T_ref = 273.15; // Temperatura de referência [K]
const sutherland_C = 110.4;      // Constante de Sutherland para o ar [K]

// Condições do ambiente
let gravity = 9.81; // Aceleração da gravidade [m/s²]
let windSpeed = 0;  // Velocidade do vento [m/s]
const ambientPressure = 101325.0; // Pressão atmosférica padrão [Pa]
const ambientTemperature = 293.15; // Temperatura ambiente [K] (~20°C)
const ambientDensity = ambientPressure / (R_specific_air * ambientTemperature); // Densidade ambiente (ρ₀) [kg/m³] - Referência para Boussinesq

let pressureStats = { min: 0, max: 0, avg: 0 };
let temperatureStats = { min: 0, max: 0, avg: 0 };
resetFields();

// --- INÍCIO DA REFATORAÇÃO: Nova função de formatação de números ---
/**
 * Formata um número para exibição na UI.
 * Usa notação de ponto fixo para valores em uma faixa "razoável".
 * Usa notação científica com HTML (ex: 1.23 &times; 10<sup>-4</sup>) para outros valores.
 * @param {number} num O número a ser formatado.
 * @param {number} precision O número de casas decimais para o coeficiente.
 * @returns {string} O número formatado como uma string HTML.
 */
function formatNumberForDisplay(num, precision = 2) {
    if (num === null || num === undefined || !isFinite(num)) {
        return "N/A";
    }
    if (Math.abs(num) < 1e-9) return "0.00";

    // Para números em uma faixa "razoável", use notação de ponto fixo.
    const absNum = Math.abs(num);
    if (absNum >= 0.001 && absNum < 10000) {
        if (absNum >= 1000) return num.toFixed(0);
        if (absNum >= 100) return num.toFixed(1);
        if (absNum >= 1) return num.toFixed(2);
        return num.toFixed(3); // Para números pequenos como 0.123
    }

    // Caso contrário, use a notação científica personalizada com HTML.
    const parts = num.toExponential(precision).split('e');
    const coefficient = parts[0];
    const exponent = parts[1].replace('+', ''); // Remove o '+' do expoente
    return `${coefficient} &times; 10<sup>${exponent}</sup>`;
}
// --- FIM DA REFATORAÇÃO ---


// --- FUNÇÕES DE SETUP ---
function resetFields() {
  u.fill(0); v.fill(0); u_prev.fill(0); v_prev.fill(0);
  r.fill(0); g.fill(0); b.fill(0); r_prev.fill(0); g_prev.fill(0); b_prev.fill(0);
  div.fill(0); p.fill(0); vort.fill(0); obstacle.fill(0);
  
  temp.fill(ambientTemperature); temp_prev.fill(ambientTemperature);
  
  // Atualiza os campos derivados (densidade, viscosidade, etc.)
  updatePhysicsFields();

  obstacles = [{
    shape: 'circle',
    x: Math.floor(NX * 0.5),
    y: Math.floor(NY * 0.75),
    radius: Math.floor(Math.min(NX, NY) * 0.12)
  }];
  updateObstacles();
}

function updateObstacles() {
    obstacle.fill(0);
    obstacles.forEach(obs => {
        const cx = obs.x, cy = obs.y, rad = obs.radius;
        const i_min = Math.max(1, Math.floor(cx - rad));
        const i_max = Math.min(NX, Math.ceil(cx + rad));
        const j_min = Math.max(1, Math.floor(cy - rad));
        const j_max = Math.min(NY, Math.ceil(cy + rad));
        
        for (let j = j_min; j <= j_max; j++) {
            for (let i = i_min; i <= i_max; i++) {
                let isInside = false;
                if (obs.shape === 'circle') {
                    const dx_grid = i - cx, dy_grid = j - cy;
                    if (dx_grid * dx_grid + dy_grid * dy_grid <= rad * rad) isInside = true;
                } else if (obs.shape === 'square') {
                    if (Math.abs(i - cx) <= rad && Math.abs(j - cy) <= rad) isInside = true;
                }
                if(isInside) obstacle[IX(i, j)] = 1;
            }
        }
    });
}

// --- SOLVER DE FÍSICA ---
function setBoundaries(field) {
    // Condições de contorno abertas (Neumann)
    for (let i = 1; i <= NX; i++) {
        field[IX(i, 0)] = field[IX(i, 1)];
        field[IX(i, NY + 1)] = field[IX(i, NY)];
    }
    for (let j = 1; j <= NY; j++) {
        field[IX(0, j)] = field[IX(1, j)];
        field[IX(NX + 1, j)] = field[IX(NX, j)];
    }
    // Cantos
    field[IX(0, 0)] = 0.5 * (field[IX(1, 0)] + field[IX(0, 1)]);
    field[IX(0, NY + 1)] = 0.5 * (field[IX(1, NY + 1)] + field[IX(0, NY)]);
    field[IX(NX + 1, 0)] = 0.5 * (field[IX(NX, 0)] + field[IX(NX + 1, 1)]);
    field[IX(NX + 1, NY + 1)] = 0.5 * (field[IX(NX, NY + 1)] + field[IX(NX + 1, NY)]);

    // Forçar velocidade zero dentro dos obstáculos (condição de não-deslizamento/impermeabilidade)
    if (field === u || field === v) {
        for (let j = 1; j <= NY; j++) {
            for (let i = 1; i <= NX; i++) {
                if (obstacle[IX(i, j)]) {
                    u[IX(i, j)] = 0;
                    v[IX(i, j)] = 0;
                }
            }
        }
    }
}


const tubeTop = Math.floor(NY * 0.42);
const tubeBottom = Math.floor(NY * 0.58);
function addInjection() {
  const x = 2;
  const injectionRatio = Math.min(1.0, (tubeBottom - tubeTop) / NY);
  for (let j = tubeTop; j <= tubeBottom; j++) {
    const idx = IX(x, j);
    if(obstacle[idx]) continue;
    
    // Adiciona cor (visual)
    r[idx] = injectionColor.r / 255;
    g[idx] = injectionColor.g / 255;
    b[idx] = injectionColor.b / 255;
    
    // Define velocidade de injeção
    u[idx] = injectionVelocity;
    
    // Define temperatura de injeção
    temp[idx] = injectionTemp;
  }
}

// --- CORREÇÃO: Função central para atualizar todos os campos físicos dependentes ---
function updatePhysicsFields() {
    for (let i = 0; i < size; i++) {
        const T = Math.max(1, temp[i]); // Evita temperatura zero ou negativa

        // 1. Equação de Estado (Gás Ideal) para calcular a densidade LOCAL ρ(T)
        // Usada *apenas* para o termo de empuxo na Aproximação de Boussinesq.
        density[i] = ambientPressure / (R_specific_air * T);

        // 2. Lei de Sutherland para viscosidade dinâmica (μ)
        visc_field[i] = base_dynamic_viscosity * Math.pow(T / sutherland_T_ref, 1.5) * (sutherland_T_ref + sutherland_C) / (T + sutherland_C);

        // 3. Difusividade Térmica (α)
        // α = k / (ρ₀ * cp). Usa-se a densidade ambiente constante (ρ₀) para consistência.
        thermal_diff_field[i] = thermalConductivity_k / (ambientDensity * specificHeat_cp);
    }
}

function applyGlobalForces(dt) {
    const windX = windSpeed;
    for (let j = 1; j <= NY; j++) {
        for (let i = 1; i <= NX; i++) {
            const idx = IX(i, j);
            if (obstacle[idx]) continue;

            // Vento
            u[idx] += dt * windX;
            
            // CORREÇÃO: Empuxo de Boussinesq
            // A força de empuxo é F = g * (ρ₀ - ρ_local)
            // A aceleração é a = F / m = F / ρ₀ = g * (ρ₀ - ρ_local) / ρ₀
            const currentDensity = density[idx];
            if (gravity > 0) {
                const buoyancy_accel = gravity * (ambientDensity - currentDensity) / ambientDensity;
                v[idx] += dt * buoyancy_accel;
            }
        }
    }
}

function applyVorticityConfinement(u_field, v_field, dt_val) {
  // NOTA: Confinamento de vorticidade é um termo artificial para adicionar detalhes visuais
  // que são perdidos pela dissipação numérica. Não é um termo físico rigoroso das
  // equações de Navier-Stokes. A força é ajustada por um parâmetro (vorticityStrength).

  for (let j = 1; j <= NY; j++) {
    for (let i = 1; i <= NX; i++) {
      const du_dy = (u_field[IX(i, j + 1)] - u_field[IX(i, j - 1)]) / (2 * dy);
      const dv_dx = (v_field[IX(i + 1, j)] - v_field[IX(i - 1, j)]) / (2 * dx);
      vort[IX(i, j)] = dv_dx - du_dy;
    }
  }

  for (let j = 2; j < NY; j++) {
    for (let i = 2; i < NX; i++) {
      const idx = IX(i, j);
      if(obstacle[idx]) continue;
      
      const gradX = (Math.abs(vort[IX(i + 1, j)]) - Math.abs(vort[IX(i - 1, j)])) / (2 * dx);
      const gradY = (Math.abs(vort[IX(i, j + 1)]) - Math.abs(vort[IX(i, j - 1)])) / (2 * dy);
      const len = Math.sqrt(gradX * gradX + gradY * gradY) + 1e-9;
      
      const force_x = (gradY / len) * vort[idx];
      const force_y = (-gradX / len) * vort[idx];
      
      // A força de confinamento é adicionada como uma aceleração.
      u_field[idx] += vorticityStrength * force_x * dt_val;
      v_field[idx] += vorticityStrength * force_y * dt_val;
    }
  }
}


// Advecção baseada nas dimensões físicas (dx, dy)
// NOTA: Este método semi-Lagrangiano é incondicionalmente estável mas não é
// formalmente conservativo. Para simulações científicas, métodos de volume finito
// (ex: MacCormack, WENO) são preferíveis.
function advect(field, field0, u_field, v_field, dt) {
  const dtx = dt / dx;
  const dty = dt / dy;

  for (let j = 1; j <= NY; j++) {
    for (let i = 1; i <= NX; i++) {
        if(obstacle[IX(i,j)]) continue;
      let x = i - dtx * u_field[IX(i, j)];
      let y = j - dty * v_field[IX(i, j)];
      
      x = Math.max(0.5, Math.min(NX + 0.5, x));
      y = Math.max(0.5, Math.min(NY + 0.5, y));

      const i0 = Math.floor(x), i1 = i0 + 1;
      const j0 = Math.floor(y), j1 = j0 + 1;
      
      const s1 = x - i0, s0 = 1 - s1;
      const t1 = y - j0, t0 = 1 - t1;

      field[IX(i, j)] = s0 * (t0 * field0[IX(i0, j0)] + t1 * field0[IX(i0, j1)]) +
                        s1 * (t0 * field0[IX(i1, j0)] + t1 * field0[IX(i1, j1)]);
    }
  }
   setBoundaries(field);
}

// CORREÇÃO: Difusão conservativa para coeficientes variáveis.
// Resolve: (I - dt * ∇·(α∇))x = x₀, onde α é o campo de difusividade.
function diffuse(field, field0, diff_field, dt) {
    const dx2 = dx * dx;
    const dy2 = dy * dy;
    
    // Aumentado para 40 iterações para melhor convergência. Para problemas mais difíceis,
    // um solver mais avançado (Multigrid, Conjugate Gradient) seria necessário.
    for (let k = 0; k < 40; k++) {
        for (let j = 1; j <= NY; j++) {
            for (let i = 1; i <= NX; i++) {
                if (obstacle[IX(i, j)]) continue;
                
                const idx = IX(i, j);
                
                // Média harmônica/aritmética dos coeficientes nas faces da célula
                const alpha_e = 0.5 * (diff_field[idx] + diff_field[IX(i + 1, j)]);
                const alpha_w = 0.5 * (diff_field[idx] + diff_field[IX(i - 1, j)]);
                const alpha_n = 0.5 * (diff_field[idx] + diff_field[IX(i, j + 1)]);
                const alpha_s = 0.5 * (diff_field[idx] + diff_field[IX(i, j - 1)]);

                const denominator = 1 + dt * (
                    (alpha_e + alpha_w) / dx2 +
                    (alpha_n + alpha_s) / dy2
                );
                if (denominator < 1e-9) continue;

                const numerator = field0[idx] + dt * (
                    (alpha_e * field[IX(i + 1, j)] + alpha_w * field[IX(i - 1, j)]) / dx2 +
                    (alpha_n * field[IX(i, j + 1)] + alpha_s * field[IX(i, j - 1)]) / dy2
                );
                
                field[idx] = numerator / denominator;
            }
        }
        setBoundaries(field);
    }
}


// CORREÇÃO: Projeção com unidades físicas corretas para garantir campo de velocidade livre de divergência.
// 1. Resolve a Equação de Poisson para a pressão: ∇²p = (ρ₀/dt) * ∇·u
// 2. Corrige a velocidade: u' = u - (dt/ρ₀) * ∇p
function project(u_field, v_field, p_field, div_field, dt) {
    // NOTA: Esta implementação em uma grade colocalizada pode sofrer de problemas de
    // acoplamento pressão-velocidade ("checkerboarding"). Uma grade escalonada (MAC)
    // é a solução padrão para evitar isso em simulações de alta fidelidade.

    // 1. Calcula a divergência ∇·u do campo de velocidade
    for (let j = 1; j <= NY; j++) {
        for (let i = 1; i <= NX; i++) {
            if (obstacle[IX(i,j)]) {
                div_field[IX(i,j)] = 0;
                continue;
            }
            div_field[IX(i, j)] = (u_field[IX(i + 1, j)] - u_field[IX(i - 1, j)]) / (2 * dx) +
                                 (v_field[IX(i, j + 1)] - v_field[IX(i, j - 1)]) / (2 * dy);
            p_field[IX(i, j)] = 0;
        }
    }
    setBoundaries(div_field);
    setBoundaries(p_field);

    // 2. Resolve a Equação de Poisson para a pressão (∇²p = (ρ₀/dt) * ∇·u) usando Gauss-Seidel
    const dx2 = dx * dx;
    const dy2 = dy * dy;
    const rho_by_dt = ambientDensity / dt;
    
    for (let k = 0; k < 80; k++) {
        for (let j = 1; j <= NY; j++) {
            for (let i = 1; i <= NX; i++) {
                if (obstacle[IX(i, j)]) continue;
                
                const numerator = (
                    (p_field[IX(i - 1, j)] + p_field[IX(i + 1, j)]) / dx2 +
                    (p_field[IX(i, j - 1)] + p_field[IX(i, j + 1)]) / dy2
                ) - (rho_by_dt * div_field[IX(i, j)]);
                
                const denominator = 2/dx2 + 2/dy2;
                
                p_field[IX(i, j)] = numerator / denominator;
            }
        }
        setBoundaries(p_field);
    }
    
    // 3. Corrige a velocidade subtraindo o gradiente de pressão (u' = u - (dt/ρ₀) * ∇p)
    const dt_by_rho = dt / ambientDensity;
    for (let j = 1; j <= NY; j++) {
        for (let i = 1; i <= NX; i++) {
            if (obstacle[IX(i, j)]) continue;
            
            const grad_p_x = (p_field[IX(i + 1, j)] - p_field[IX(i - 1, j)]) / (2 * dx);
            const grad_p_y = (p_field[IX(i, j + 1)] - p_field[IX(i, j - 1)]) / (2 * dy);
            
            u_field[IX(i, j)] -= dt_by_rho * grad_p_x;
            v_field[IX(i, j)] -= dt_by_rho * grad_p_y;
        }
    }
    setBoundaries(u_field);
    setBoundaries(v_field);
}

function calculateStats(field, statsObject) {
  let min = Infinity, max = -Infinity, sum = 0, count = 0;
  for (let j = 1; j <= NY; j++) {
    for (let i = 1; i <= NX; i++) {
      if (!obstacle[IX(i, j)]) {
        const value = field[IX(i,j)];
        min = Math.min(min, value); max = Math.max(max, value);
        sum += value; count++;
      }
    }
  }
  statsObject.min = min; statsObject.max = max; statsObject.avg = count > 0 ? sum / count : 0;
}

// --- NOVOS CÁLCULOS FÍSICOS ---
function calculateKineticEnergyHistogram() {
    const bins = 20;
    const histogram = new Array(bins).fill(0);
    const cellVolume = dx * dy; // Volume da célula (assumindo profundidade de 1m)
    
    let maxKE = 0;
    const kineticEnergies = [];
    
    // Primeiro, calcular todas as energias cinéticas
    for (let j = 1; j <= NY; j++) {
        for (let i = 1; i <= NX; i++) {
            const idx = IX(i, j);
            if (!obstacle[idx]) {
                const speed2 = u[idx] * u[idx] + v[idx] * v[idx];
                const ke = 0.5 * ambientDensity * cellVolume * speed2; // E = 0.5 * m * v² (usa ρ₀)
                kineticEnergies.push(ke);
                maxKE = Math.max(maxKE, ke);
            }
        }
    }
    
    // Criar histograma
    const binSize = maxKE > 1e-9 ? maxKE / bins : 1;
    kineticEnergies.forEach(ke => {
        const binIndex = Math.min(Math.floor(ke / binSize), bins - 1);
        histogram[binIndex]++;
    });
    
    return {
        histogram: histogram,
        binSize: binSize,
        maxKE: maxKE,
        totalCells: kineticEnergies.length
    };
}

function calculateTotalThermalEnergy() {
    let totalThermalEnergy = 0;
    const cellVolume = dx * dy; // Volume da célula (assumindo profundidade de 1m)
    
    for (let j = 1; j <= NY; j++) {
        for (let i = 1; i <= NX; i++) {
            const idx = IX(i, j);
            if (!obstacle[idx]) {
                // Energia térmica = m * cp * (T - T_ref)
                const mass = density[idx] * cellVolume; // Usa densidade local para massa de calor
                const thermalEnergy = mass * specificHeat_cp * (temp[idx] - ambientTemperature);
                totalThermalEnergy += thermalEnergy;
            }
        }
    }
    
    return totalThermalEnergy;
}

function calculateTotalMomentum() {
    let momentumX = 0;
    let momentumY = 0;
    const cellVolume = dx * dy; // Volume da célula (assumindo profundidade de 1m)
    
    for (let j = 1; j <= NY; j++) {
        for (let i = 1; i <= NX; i++) {
            const idx = IX(i, j);
            if (!obstacle[idx]) {
                const mass = ambientDensity * cellVolume; // Momento usa ρ₀
                momentumX += mass * u[idx]; // p = m * v
                momentumY += mass * v[idx];
            }
        }
    }
    
    return { x: momentumX, y: momentumY };
}


// --- PASSO DA SIMULAÇÃO (REORGANIZADO PARA CLAREZA FÍSICA) ---
let dt; // dt é definido no loop de simulação
function step() {
  dt = parseFloat(document.getElementById('dtRange').value) * 0.005;

  // 1. Atualiza propriedades derivadas
  updatePhysicsFields();

  // 2. Aplica forças externas (gravidade/empuxo Boussinesq, vento)
  applyGlobalForces(dt);

  // 3. Adiciona fontes (injeção de fluido)
  if (injectionOn) addInjection();

  // ------------ VELOCIDADES ------------
  // Antes de advect, salve o estado atual em u_prev/v_prev (fonte para advect)
  u_prev.set(u);
  v_prev.set(v);

  // 4. Advecção do campo de velocidade (transporta quantidade de movimento)
  advect(u, u_prev, u_prev, v_prev, dt);
  advect(v, v_prev, u_prev, v_prev, dt);

  // 5. Aplica confinamento de vorticidade
  if (vorticityStrength > 0) applyVorticityConfinement(u, v, dt);

  // 6. Difusão da velocidade (usar kinematic_visc_field)
  // Salvar estado atual antes da difusão (fonte)
  u_prev.set(u);
  v_prev.set(v);

  const kinematic_visc_field = new Float32Array(size);
  for (let i = 0; i < size; i++) {
    kinematic_visc_field[i] = visc_field[i] / ambientDensity;
  }
  diffuse(u, u_prev, kinematic_visc_field, dt);
  diffuse(v, v_prev, kinematic_visc_field, dt);

  // 7. Projeção para garantir incompressibilidade (∇·u = 0)
  project(u, v, p, div, dt);

  // ------------ TEMPERATURA ------------
  // Salvar temp atual, advect/diffuse sobre essa fonte
  temp_prev.set(temp);
  advect(temp, temp_prev, u, v, dt);
  // Para difusão térmica, salve novamente (fonte = resultado da advect)
  temp_prev.set(temp);
  diffuse(temp, temp_prev, thermal_diff_field, dt);

  // ------------ CORES (VISUAL) ------------
  r_prev.set(r); g_prev.set(g); b_prev.set(b);
  advect(r, r_prev, u, v, dt);
  advect(g, g_prev, u, v, dt);
  advect(b, b_prev, u, v, dt);

  // 11. Estatísticas e gravação (inalterado)
  calculateStats(p, pressureStats);
  calculateStats(temp, temperatureStats);

  if (isRecording) {
    const prevSimTime = simulationTimeRecorded;
    simulationTimeRecorded += dt;
    const sampleInterval = 1.0 / recordingFrequency;
    const prevSamplePoint = Math.floor(prevSimTime / sampleInterval);
    const currentSamplePoint = Math.floor(simulationTimeRecorded / sampleInterval);
    if (currentSamplePoint > prevSamplePoint || originalRecordedData.length === 0) {
      collectFrameData(dt);
    }
    if (simulationTimeRecorded >= recordingDuration) {
      stopRecording();
    }
  }
}


// --- VISUALIZAÇÃO ---
let visualizationMode = 'density';
let showVel = false, showTrails = true;

function getViridisColor(value, min, max) {
  const t = Math.max(0, Math.min(1, (value - min) / (max - min)));

  // 6 pontos da escala Viridis
  const viridis = [
    [68, 1, 84],     // roxo escuro
    [59, 82, 139],   // azul
    [33, 145, 140],  // verde água
    [94, 201, 98],   // verde
    [253, 231, 37]   // amarelo
  ];

  const idx = t * (viridis.length - 1);
  const i0 = Math.floor(idx);
  const i1 = Math.min(i0 + 1, viridis.length - 1);
  const f = idx - i0;

  const r = Math.round(viridis[i0][0] + f * (viridis[i1][0] - viridis[i0][0]));
  const g = Math.round(viridis[i0][1] + f * (viridis[i1][1] - viridis[i0][1]));
  const b = Math.round(viridis[i0][2] + f * (viridis[i1][2] - viridis[i0][2]));

  return { r, g, b };
}


function getDivergingColor(value, scale) {
  const t = Math.max(-1, Math.min(1, value / scale));
  if (t < 0) {
    const intensity = Math.abs(t);
    return { r: Math.round(255 * intensity), g: Math.round(100 * intensity), b: Math.round(255 * intensity) };
  } else {
    return { r: Math.round(255 * t), g: Math.round(165 * t), b: 0 };
  }
}

function render() {
  const img = ctx.createImageData(canvas.width, canvas.height);
  const data = img.data;
  const sx = canvas.width / NX, sy = canvas.height / NY;

  if (!showTrails) data.fill(0);

  for (let j = 1; j <= NY; j++) {
    for (let i = 1; i <= NX; i++) {
      const idx = IX(i, j);
      const x0 = Math.floor((i - 1) * sx), y0 = Math.floor((j - 1) * sy);
      const x1 = Math.floor(i * sx), y1 = Math.floor(j * sy);
      
      let R = 0, G = 0, B = 0;
      if (obstacle[idx]) {
        R = G = B = 80;
      } else {
        if (visualizationMode === 'density') {
          R = Math.round(255 * r[idx]); G = Math.round(255 * g[idx]); B = Math.round(255 * b[idx]);
        }
        else if (visualizationMode === 'temperature') {
            const range = (temperatureStats.max - temperatureStats.min) / temperatureSensitivity;
            const mid = (temperatureStats.max + temperatureStats.min) / 2;
            const minScaled = mid - range / 2;
            const maxScaled = mid + range / 2;
            const c = getViridisColor(temp[idx], minScaled, maxScaled);
            R = c.r; G = c.g; B = c.b;
        }
        else if (visualizationMode === 'pressure') {
            const centered = (p[idx] - pressureStats.avg) * pressureSensitivity + pressureStats.avg;
            const c = getViridisColor(centered, pressureStats.min, pressureStats.max);
            R = c.r; G = c.g; B = c.b;
        }
        else if (visualizationMode === 'vorticity') { 
            const c = getDivergingColor(vort[idx], vorticityVisScale);
            R=c.r; G=c.g; B=c.b; 
        }
      }
      for (let y = y0; y < y1; y++) {
        for (let x = x0; x < x1; x++) {
          const pIdx = (y * canvas.width + x) * 4;
          data[pIdx] = R; data[pIdx + 1] = G; data[pIdx + 2] = B; data[pIdx + 3] = 255;
        }
      }
    }
  }
  ctx.putImageData(img, 0, 0);

  if (injectionOn) {
    ctx.fillStyle = 'rgba(255, 80, 0, 0.6)';
    ctx.fillRect(0, tubeTop * sy, 10, (tubeBottom - tubeTop + 1) * sy);
  }

  if (showPressureContours && visualizationMode === 'pressure') drawPressureContours(pressureStats.min, pressureStats.max);
  if (showVel) drawVelocityVectors(sx, sy);
}

function drawVelocityVectors(sx, sy) {
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)'; ctx.lineWidth = 1;
    const step = 8;
    for (let j = step; j <= NY; j += step) {
      for (let i = step; i <= NX; i += step) {
        const idx = IX(i, j);
        const speed = Math.hypot(u[idx], v[idx]);
        if (obstacle[idx] || speed < 0.01) continue;
        const x = (i - 0.5) * sx; 
        const y = (j - 0.5) * sy;
        const scale = 2.0; // Fator de escala para visualização
        ctx.beginPath(); 
        ctx.moveTo(x, y); 
        ctx.lineTo(x + u[idx] * scale, y + v[idx] * scale); 
        ctx.stroke();
      }
    }
}

function drawPressureContours(minP, maxP) {
  ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)'; // Cor sutil para os contornos
  ctx.lineWidth = 1;
  
  const levels = 6; // Número de linhas de contorno a serem desenhadas
  const sx = canvas.width / NX;
  const sy = canvas.height / NY;
  
  for (let level = 1; level < levels; level++) {
    const threshold = minP + (maxP - minP) * level / levels;
    
    for (let j = 0; j < NY; j++) {
      for (let i = 0; i < NX; i++) {
        if (obstacle[IX(i, j)]) continue;
        
        const corners = [ p[IX(i, j)], p[IX(i + 1, j)], p[IX(i + 1, j + 1)], p[IX(i, j + 1)] ];
        let caseIndex = 0;
        if (corners[0] > threshold) caseIndex |= 1; 
        if (corners[1] > threshold) caseIndex |= 2;
        if (corners[2] > threshold) caseIndex |= 4; 
        if (corners[3] > threshold) caseIndex |= 8;

        const startX = i * sx, startY = j * sy;
        const midX = startX + sx * 0.5, midY = startY + sy * 0.5;

        // Interpolação linear para encontrar pontos exatos na borda da célula
        const t_top = (threshold - corners[0]) / (corners[1] - corners[0]);
        const t_bot = (threshold - corners[3]) / (corners[2] - corners[3]);
        const t_lft = (threshold - corners[0]) / (corners[3] - corners[0]);
        const t_rgt = (threshold - corners[1]) / (corners[2] - corners[1]);

        const mid = {
            top:    { x: startX + sx * t_top, y: startY },
            right:  { x: startX + sx, y: startY + sy * t_rgt },
            bottom: { x: startX + sx * t_bot, y: startY + sy },
            left:   { x: startX, y: startY + sy * t_lft }
        };
        
        ctx.beginPath();
        switch (caseIndex) {
            case 1: case 14: ctx.moveTo(mid.left.x, mid.left.y); ctx.lineTo(mid.bottom.x, mid.bottom.y); break;
            case 2: case 13: ctx.moveTo(mid.right.x, mid.right.y); ctx.lineTo(mid.bottom.x, mid.bottom.y); break;
            case 3: case 12: ctx.moveTo(mid.left.x, mid.left.y); ctx.lineTo(mid.right.x, mid.right.y); break;
            case 4: case 11: ctx.moveTo(mid.top.x, mid.top.y); ctx.lineTo(mid.right.x, mid.right.y); break;
            case 5: ctx.moveTo(mid.left.x, mid.left.y); ctx.lineTo(mid.top.x, mid.top.y); ctx.moveTo(mid.right.x, mid.right.y); ctx.lineTo(mid.bottom.x, mid.bottom.y); break;
            case 6: case 9: ctx.moveTo(mid.top.x, mid.top.y); ctx.lineTo(mid.bottom.x, mid.bottom.y); break;
            case 7: case 8: ctx.moveTo(mid.left.x, mid.left.y); ctx.lineTo(mid.top.x, mid.top.y); break;
            case 10: ctx.moveTo(mid.left.x, mid.left.y); ctx.lineTo(mid.bottom.x, mid.bottom.y); ctx.moveTo(mid.top.x, mid.top.y); ctx.lineTo(mid.right.x, mid.right.y); break;
        }
        ctx.stroke();
      }
    }
  }
}

// --- LOOP & UI ---
let lastTime = 0, frameCount = 0;
function simLoop(t) {
  const elapsed = t - lastTime;
  if (elapsed > 1000) {
      document.getElementById('perfStats').textContent = `${(frameCount * 1000 / elapsed).toFixed(1)}`;
      lastTime = t; frameCount = 0;
  }
  
  if (!paused) step();
  
  if (isRecording) {
    const clampedSimTime = Math.max(0, Math.min(simulationTimeRecorded, recordingDuration));
    const fraction = recordingDuration > 0 ? (clampedSimTime / recordingDuration) : 1;
    const progressEl = recordBtn._progressEl || document.getElementById('recordProgress');
    if (progressEl) { 
        progressEl.style.width = (fraction * 100).toFixed(2) + '%'; 
    }
    const span = recordBtn.querySelector('.recordText');
    if (span) { 
        span.textContent = `Gravando... ${clampedSimTime.toFixed(1)}s / ${recordingDuration.toFixed(1)}s`;
    }
  }


  render(); updateUI(); frameCount++;
  requestAnimationFrame(simLoop);
}

function updateUI() {
  const modeNames = {'density': 'Corante (Visual)', 'pressure': 'Pressão (Pa)', 'vorticity': 'Vorticidade', 'temperature': 'Temperatura'};
  document.getElementById('currentMode').textContent = modeNames[visualizationMode];
  document.getElementById('simState').textContent = paused ? 'Pausado' : 'Executando';
  document.getElementById('statusIndicator').className = `status-indicator ${paused ? 'status-paused' : 'status-running'}`;
  
  const pressureInfo = document.getElementById('pressureInfo');
  const temperatureInfo = document.getElementById('temperatureInfo');
  const colorScale = document.getElementById('colorScale');
  
  pressureInfo.style.display = 'none';
  temperatureInfo.style.display = 'none';
  colorScale.style.display = 'none';

  if (visualizationMode === 'pressure') {
    pressureInfo.style.display = 'block';
    colorScale.style.display = 'block';
    document.getElementById('scaleTitle').textContent = "Pressão (Pa)";
    // --- INÍCIO DA REFATORAÇÃO: Usar nova função de formatação e innerHTML ---
    document.getElementById('maxPressure').innerHTML = formatNumberForDisplay(pressureStats.max, 2);
    document.getElementById('minPressure').innerHTML = formatNumberForDisplay(pressureStats.min, 2);
    document.getElementById('avgPressure').innerHTML = formatNumberForDisplay(pressureStats.avg, 2);
    document.getElementById('minScale').innerHTML = formatNumberForDisplay(pressureStats.min, 1);
    document.getElementById('maxScale').innerHTML = formatNumberForDisplay(pressureStats.max, 1);
    // --- FIM DA REFATORAÇÃO ---
  } else if (visualizationMode === 'temperature') {
    temperatureInfo.style.display = 'block';
    colorScale.style.display = 'block';
    document.getElementById('scaleTitle').textContent = "Escala de Temperatura (K)";
    // --- INÍCIO DA REFATORAÇÃO: Usar nova função de formatação ---
    document.getElementById('maxTemp').innerHTML = formatNumberForDisplay(temperatureStats.max, 2);
    document.getElementById('minTemp').innerHTML = formatNumberForDisplay(temperatureStats.min, 2);
    document.getElementById('avgTemp').innerHTML = formatNumberForDisplay(temperatureStats.avg, 2);
    document.getElementById('minScale').innerHTML = formatNumberForDisplay(temperatureStats.min, 2);
    document.getElementById('maxScale').innerHTML = formatNumberForDisplay(temperatureStats.max, 2);
    // --- FIM DA REFATORAÇÃO ---
  }
}

requestAnimationFrame(simLoop);

// --- EVENT LISTENERS ---
function setupSlider(rangeId, valId, callback, formatter = v => v) {
    const range = document.getElementById(rangeId);
    const val = document.getElementById(valId);
    if (!range || !val) {
        console.warn(`Elemento de slider não encontrado: ${rangeId} ou ${valId}`);
        return;
    }
    const update = () => {
        // --- INÍCIO DA REFATORAÇÃO: Usar innerHTML para permitir formatação rica ---
        val.innerHTML = formatter(range.value);
        // --- FIM DA REFATORAÇÃO ---
        if (callback) callback(parseFloat(range.value));
    };
    range.addEventListener('input', update);
    update();
}

setupSlider('dtRange', 'dtVal');
setupSlider('injectRange', 'injectVal', v => injectionVelocity = v);
setupSlider('injectTempRange', 'injectTempVal', v => injectionTemp = v);
// --- INÍCIO DA REFATORAÇÃO: Atualizar o formatador para a nova função ---
setupSlider('viscRange', 'viscVal', v => base_dynamic_viscosity = v, v => formatNumberForDisplay(parseFloat(v), 1));
// --- FIM DA REFATORAÇÃO ---
setupSlider('vortRange', 'vortVal', v => vorticityStrength = v);
setupSlider('vortSensRange', 'vortSensVal', v => vorticityVisScale = v);
setupSlider('pressureSensRange', 'pressureSensVal', v => pressureSensitivity = v);
setupSlider('gravityRange', 'gravityVal', v => gravity = v);
setupSlider('windSpeedRange', 'windSpeedVal', v => windSpeed = v);

document.getElementById('velToggle').addEventListener('change', e => showVel = e.target.checked);
document.getElementById('trailToggle').addEventListener('change', e => showTrails = e.target.checked);
document.getElementById('pressureContours').addEventListener('change', e => showPressureContours = e.target.checked);
document.getElementById('visModeSelect').addEventListener('change', e => visualizationMode = e.target.value);
document.getElementById('addCircleBtn').addEventListener('click', () => {
    obstacles.push({ shape: 'circle', x: Math.floor(NX * 0.25), y: Math.floor(NY * 0.5), radius: Math.floor(Math.min(NX, NY) * 0.08) });
    updateObstacles();
});
document.getElementById('addSquareBtn').addEventListener('click', () => {
    obstacles.push({ shape: 'square', x: Math.floor(NX * 0.75), y: Math.floor(NY * 0.5), radius: Math.floor(Math.min(NX, NY) * 0.08) });
    updateObstacles();
});
document.getElementById('clearObstaclesBtn').addEventListener('click', () => {
    obstacles = [];
    updateObstacles();
});

document.getElementById('pauseBtn').addEventListener('click', () => {
  paused = !paused;
  if (!paused) {
    canvasTooltip.style.display = 'none';
  }
  const btn = document.getElementById('pauseBtn');
  btn.innerHTML = paused ? 
    `<svg class="icon" viewBox="0 0 24 24"><path d="M8 5v14l11-7z"/></svg> Retomar` : 
    `<svg class="icon" viewBox="0 0 24 24"><path d="M6 19h4V5H6v14zm8-14v14h4V5h-4z"/></svg> Pausar`;
});
document.getElementById('stepBtn').addEventListener('click', () => { if (paused) step(); });
document.getElementById('resetBtn').addEventListener('click', resetFields);
document.getElementById('toggleWind').addEventListener('click', () => injectionOn = !injectionOn);

function hexToRgb(hex) {
    const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? { r: parseInt(result[1], 16), g: parseInt(result[2], 16), b: parseInt(result[3], 16) } : null;
}
document.getElementById('colorPicker').addEventListener('input', e => injectionColor = hexToRgb(e.target.value));

// Interação do mouse
let dragging = false, dragMode = null;
function getMouseCell(ev) {
    const rect = canvas.getBoundingClientRect();
    const i = Math.floor((ev.clientX - rect.left) / (canvas.width / NX)) + 1;
    const j = Math.floor((ev.clientY - rect.top) / (canvas.height / NY)) + 1;
    if (i >= 1 && i <= NX && j >= 1 && j <= NY) return { i, j };
    return null;
}

function getObstacleAt(i, j) {
    for (let k = obstacles.length - 1; k >= 0; k--) {
        const obs = obstacles[k];
        const dx_grid = i - obs.x;
        const dy_grid = j - obs.y;
        if (obs.shape === 'circle') {
            if (dx_grid * dx_grid + dy_grid * dy_grid <= obs.radius * obs.radius) return k;
        } else if (obs.shape === 'square') {
            if (Math.abs(dx_grid) <= obs.radius && Math.abs(dy_grid) <= obs.radius) return k;
        }
    }
    return -1;
}

canvas.addEventListener('mousedown', (ev) => {
  const cell = getMouseCell(ev);
  if (!cell) return;
  dragging = true;
  draggedObstacleIndex = getObstacleAt(cell.i, cell.j);
  
  if (ev.shiftKey) {
      dragMode = 'gust'; 
  } else if (ev.ctrlKey) {
      dragMode = 'drain'; 
  } else if (ev.button === 2) { 
      dragMode = 'vortex'; 
  } else if (draggedObstacleIndex !== -1) { 
      dragMode = 'obstacle'; 
  } else { 
      dragMode = 'puff'; 
  }
});

canvas.addEventListener('mousemove', (ev) => {
  if (paused) {
    const cell = getMouseCell(ev);
    if (cell) {
        const idx = IX(cell.i, cell.j);
        let tooltipContent = '';
        if (obstacle[idx]) {
            tooltipContent = 'Obstáculo';
        } else {
            switch (visualizationMode) {
                case 'density':
                    const r_val = Math.round(r[idx] * 255);
                    const g_val = Math.round(g[idx] * 255);
                    const b_val = Math.round(b[idx] * 255);
                    tooltipContent = `Cor: <b class="monospace-font">RGB(${r_val}, ${g_val}, ${b_val})</b>`;
                    break;
                case 'temperature':
                    // --- INÍCIO DA REFATORAÇÃO: Usar nova função de formatação ---
                    tooltipContent = `Temp: <b class="monospace-font">${formatNumberForDisplay(temp[idx], 2)} K</b>`;
                    // --- FIM DA REFATORAÇÃO ---
                    break;
                case 'pressure':
                    // --- INÍCIO DA REFATORAÇÃO: Usar nova função de formatação ---
                    tooltipContent = `Pressão: <b class="monospace-font">${formatNumberForDisplay(p[idx], 3)} Pa</b>`;
                    // --- FIM DA REFATORAÇÃO ---
                    break;
                case 'vorticity':
                    // --- INÍCIO DA REFATORAÇÃO: Usar nova função de formatação ---
                    tooltipContent = `Vorticidade: <b class="monospace-font">${formatNumberForDisplay(vort[idx], 4)}</b>`;
                    // --- FIM DA REFATORAÇÃO ---
                    break;
            }
        }
        canvasTooltip.innerHTML = tooltipContent;
        canvasTooltip.style.display = 'block';
        const rect = canvas.getBoundingClientRect();
        canvasTooltip.style.left = `${ev.clientX - rect.left + 15}px`;
        canvasTooltip.style.top = `${ev.clientY - rect.top + 15}px`;
    } else {
        canvasTooltip.style.display = 'none';
    }
  }

  if (!dragging) return;
  const cell = getMouseCell(ev);
  if (!cell) return;
  
  const idx = IX(cell.i, cell.j);
  if (obstacle[idx] && dragMode !== 'obstacle') return;

  if (dragMode === 'puff') {
    const color = getViridisColor(Math.random(), 0, 1);
    r[idx] = Math.min(1, r[idx] + color.r / 255 * 0.5); 
    g[idx] = Math.min(1, g[idx] + color.g / 255 * 0.5); 
    b[idx] = Math.min(1, b[idx] + color.b / 255 * 0.5);
    u[idx] += (ev.movementX || 0) * 0.1; 
    v[idx] += (ev.movementY || 0) * 0.1;
    temp[idx] = Math.min(2000, temp[idx] + 50);
  } else if (dragMode === 'obstacle' && draggedObstacleIndex !== -1) {
    obstacles[draggedObstacleIndex].x = cell.i;
    obstacles[draggedObstacleIndex].y = cell.j;
    updateObstacles();
  } else if (dragMode === 'gust' || dragMode === 'drain') {
    ev.preventDefault();
    const radius = 25;
    const strength = 5.0;

    for (let j_offset = -radius; j_offset <= radius; j_offset++) {
      for (let i_offset = -radius; i_offset <= radius; i_offset++) {
        const dist_sq = i_offset * i_offset + j_offset * j_offset;
        
        if (dist_sq > radius * radius || dist_sq < 1) continue;

        const target_i = cell.i + i_offset;
        const target_j = cell.j + j_offset;

        if (target_i <= 0 || target_i > NX || target_j <= 0 || target_j > NY) continue;

        const target_idx = IX(target_i, target_j);
        if (obstacle[target_idx]) continue;

        const dist = Math.sqrt(dist_sq);
        
        let force_x = i_offset / dist;
        let force_y = j_offset / dist;
        
        if (dragMode === 'drain') {
            force_x *= -1;
            force_y *= -1;
        }
        
        const falloff = strength / (1 + 0.1 * dist_sq);

        if (dt) {
            u[target_idx] += force_x * falloff * dt;
            v[target_idx] += force_y * falloff * dt;
        }
      }
    }
  }
});


canvas.addEventListener('mouseup', () => {
    dragging = false;
    dragMode = null;
    draggedObstacleIndex = -1;
});
canvas.addEventListener('contextmenu', (e) => e.preventDefault());

canvas.addEventListener('mouseleave', () => {
    canvasTooltip.style.display = 'none';
});

canvas.addEventListener('wheel', (ev) => {
    const cell = getMouseCell(ev);
    if (cell) {
        const obstacleIndex = getObstacleAt(cell.i, cell.j);
        if (obstacleIndex !== -1) {
            ev.preventDefault();
            obstacles[obstacleIndex].radius = Math.max(3, obstacles[obstacleIndex].radius + (ev.deltaY < 0 ? 1 : -1));
            updateObstacles();
        }
    }
}, { passive: false });

window.addEventListener('resize', () => {
  canvasW = Math.floor(window.innerWidth > 1200 ? Math.min(window.innerWidth - 420, 1280) : window.innerWidth - 40);
  canvasH = Math.floor(window.innerHeight - 80);
  canvas.width = canvasW; canvas.height = canvasH;
  updateGridDimensions();
  alloc(); resetFields();
});

// Inicializar valores da UI
showVel = document.getElementById('velToggle').checked;
showTrails = document.getElementById('trailToggle').checked;
injectionColor = hexToRgb(document.getElementById('colorPicker').value);
updateHistoryUI();

// --- SEÇÃO DE GRAVAÇÃO E ANÁLISE ---
const recordBtn = document.getElementById('recordBtn');
const recordDurationInput = document.getElementById('recordDurationInput');
const recordFrequencyInput = document.getElementById('recordFrequencyInput');
const reportModal = document.getElementById('reportModal');
const closeReportBtn = document.getElementById('closeReportBtn');
const tableSearchInput = document.getElementById('tableSearchInput');
const pdfOptionsModal = document.getElementById('pdfOptionsModal');
const csvOptionsModal = document.getElementById('csvOptionsModal');
const downloadPdfBtn = document.getElementById('downloadPdfBtn');
const downloadCsvBtn = document.getElementById('downloadCsvBtn');


function startRecording() {
    if (isRecording) return;
    
    recordingDuration = parseFloat(recordDurationInput.value);
    recordingFrequency = parseFloat(recordFrequencyInput.value);

    if (isNaN(recordingDuration) || recordingDuration <= 0) {
        alert("Por favor, insira uma duração de gravação válida.");
        return;
    }
     if (isNaN(recordingFrequency) || recordingFrequency <= 0) {
        alert("Por favor, insira uma frequência de amostragem válida.");
        return;
    }

    simulationParametersSnapshot = {
        'Resolução da Grade': `${NX} x ${NY}`,
        'Dimensões do Domínio (m)': `${DOMAIN_WIDTH.toFixed(2)} x ${DOMAIN_HEIGHT.toFixed(2)}`,
        'Passo de Tempo Médio (s)': dt ? dt.toExponential(2) : 'N/A', // Mantido 'e' para consistência interna
        'Viscosidade Din. Base (Pa·s)': base_dynamic_viscosity.toExponential(2),
        'Velocidade de Injeção (m/s)': injectionVelocity,
        'Temperatura de Injeção (K)': injectionTemp,
        'Confinamento de Vorticidade (ε)': vorticityStrength,
    };
    
    originalRecordedData = []; 
    isRecording = true; 
    simulationTimeRecorded = 0;

    const totalSec = recordingDuration.toFixed(1);
    recordBtn.innerHTML = `<svg class="icon" viewBox="0 0 24 24" style="animation: pulse 1.5s infinite; fill: #E53935;"><path d="M12 2C17.5 2 22 6.5 22 12S17.5 22 12 22 2 17.5 2 12 6.5 2 12 2M12 4C16.4 4 20 7.6 20 12S16.4 20 12 20 4 16.4 4 12 7.6 4 12 4M12 6C15.3 6 18 8.7 18 12S15.3 18 12 18 6 15.3 6 12 8.7 6 12 6Z"></path></svg>
                           <span class="recordText" style="margin-left:8px">Gravando... 0.0s / ${totalSec}s</span>`;
    recordBtn.disabled = true; 
    recordDurationInput.disabled = true;
    recordFrequencyInput.disabled = true;
    document.getElementById('reportContainer').innerHTML = '';

    let progressEl = document.getElementById('recordProgress');
    if (progressEl) {
        progressEl.style.width = '0%';
        progressEl.style.display = 'block';
        recordBtn._progressEl = progressEl;
    }
}


function stopRecording() {
    if (!isRecording && simulationTimeRecorded === 0) return;
    isRecording = false;

    recordBtn.innerHTML = `<svg class="icon" viewBox="0 0 24 24"><path fill="currentColor" d="M12,2A10,10 0 0,0 2,12A10,10 0 0,0 12,22A10,10 0 0,0 22,12A10,10 0 0,0 12,2M12,4A8,8 0 0,1 20,12A8,8 0 0,1 12,20A8,8 0 0,1 4,12A8,8 0 0,1 12,4Z"></path></svg> Iniciar Gravação`;
    recordBtn.disabled = false; 
    recordDurationInput.disabled = false;
    recordFrequencyInput.disabled = false;

    const progressEl = recordBtn._progressEl || document.getElementById('recordProgress');
    if (progressEl) {
        progressEl.style.width = '0%';
        delete recordBtn._progressEl;
    }
    
    recordedData = [...originalRecordedData];

    if (recordedData.length > 1) {
        const historyEntry = {
            timestamp: new Date(),
            parameters: JSON.parse(JSON.stringify(simulationParametersSnapshot)),
            data: [...originalRecordedData]
        };
        recordingHistory.unshift(historyEntry);
        updateHistoryUI();

        populateReportModal();
        reportModal.style.display = 'flex';
    } else {
        document.getElementById('reportContainer').innerHTML = `<small style="color: var(--error-color);">Gravação muito curta, nenhum dado foi gerado.</small>`;
    }
    simulationTimeRecorded = 0;
}


function collectFrameData(dt_val) {
    let kineticEnergy = 0, totalVorticity = 0, cellCount = 0;
    const cellVolume = dx * dy;
    
    let max_cfl = 0;
    let divergence_l2_norm_sq = 0;
    let max_vel = 0;

    for (let j = 1; j <= NY; j++) {
        for (let i = 1; i <= NX; i++) {
            const idx = IX(i,j);
            if (!obstacle[idx]) {
                const speed_u = Math.abs(u[idx]);
                const speed_v = Math.abs(v[idx]);
                const speed_sq = speed_u * speed_u + speed_v * speed_v;
                max_vel = Math.max(max_vel, Math.sqrt(speed_sq));

                const cfl = dt_val * (speed_u / dx + speed_v / dy);
                if (cfl > max_cfl) max_cfl = cfl;

                const divergence_val = div[idx];
                divergence_l2_norm_sq += divergence_val * divergence_val;
                
                kineticEnergy += 0.5 * speed_sq * ambientDensity * cellVolume;
                totalVorticity += Math.abs(vort[idx]);
                cellCount++;
            }
        }
    }
    const divergence_l2_norm = Math.sqrt(divergence_l2_norm_sq / cellCount);
    
    const thermalEnergy = calculateTotalThermalEnergy();
    const momentum = calculateTotalMomentum();
    const keHistogram = calculateKineticEnergyHistogram();
    
    originalRecordedData.push({
        'Tempo (s)': parseFloat(simulationTimeRecorded.toFixed(3)),
        'Energia Cinética (J)': kineticEnergy,
        'Energia Térmica (J)': thermalEnergy,
        'Momento Linear X (kg·m/s)': momentum.x,
        'Momento Linear Y (kg·m/s)': momentum.y,
        'Vorticidade Total (1/s)': totalVorticity,
        'Pressão Média (Pa)': pressureStats.avg,
        'Pressão Máx (Pa)': pressureStats.max,
        'Pressão Mín (Pa)': pressureStats.min,
        'Número CFL Máx': max_cfl,
        'Velocidade Máx (m/s)': max_vel,
        'Residual Div. (L2)': divergence_l2_norm,
        'Histograma EC': keHistogram
    });
}

function updateHistoryUI() {
    const historyList = document.getElementById('recordingHistoryList');
    historyList.innerHTML = '';

    if (recordingHistory.length === 0) {
        historyList.innerHTML = `<small>Nenhuma gravação no histórico.</small>`;
        return;
    }

    recordingHistory.forEach((entry, index) => {
        const li = document.createElement('li');
        li.dataset.index = index;
        
        const date = entry.timestamp.toLocaleDateString('pt-BR');
        const time = entry.timestamp.toLocaleTimeString('pt-BR', { hour: '2-digit', minute: '2-digit' });
        
        li.innerHTML = `
            <div class="history-item-content">
                <span>Gravação de ${date}</span>
                <small class="monospace-font">${time} - ${entry.data.length} amostras</small>
            </div>
        `;

        li.addEventListener('click', (e) => {
            const historyIndex = parseInt(e.currentTarget.dataset.index);
            const selectedEntry = recordingHistory[historyIndex];

            if (selectedEntry) {
                recordedData = selectedEntry.data;
                originalRecordedData = selectedEntry.data;
                simulationParametersSnapshot = selectedEntry.parameters;

                populateReportModal();
                reportModal.style.display = 'flex';
            }
        });

        historyList.appendChild(li);
    });
}


function getAggregateStats(data, key) {
    if (!data || data.length === 0 || !data[0].hasOwnProperty(key)) return { min: 0, max: 0, avg: 0, median: 0, stdDev: 0, integral: 0 };
    const values = data.map(d => d[key]); const sum = values.reduce((a, b) => a + b, 0); const avg = sum / values.length;
    const sortedValues = [...values].sort((a, b) => a - b); const mid = Math.floor(sortedValues.length / 2);
    const median = sortedValues.length % 2 !== 0 ? sortedValues[mid] : (sortedValues[mid - 1] + sortedValues[mid]) / 2;
    const variance = values.reduce((acc, val) => acc + Math.pow(val - avg, 2), 0) / values.length; const stdDev = Math.sqrt(variance);
    let integral = 0;
    for (let i = 1; i < data.length; i++) {
        const dt_report = data[i]['Tempo (s)'] - data[i-1]['Tempo (s)'];
        integral += (data[i][key] + data[i-1][key]) / 2 * dt_report;
    }
    return { min: Math.min(...values), max: Math.max(...values), avg: avg, median: median, stdDev: stdDev, integral: integral };
}

const METRIC_CONFIG = {
    'Energia Cinética (J)': { label: 'Energia Cinética Total', color: '#007AFF', icon: 'fas fa-bolt' },
    'Energia Térmica (J)': { label: 'Energia Térmica Total', color: '#e84306ff', icon: 'fas fa-thermometer-half' },
    'Momento Linear X (kg·m/s)': { label: 'Momento Linear X', color: '#28A745', icon: 'fas fa-arrow-right' },
    'Momento Linear Y (kg·m/s)': { label: 'Momento Linear Y', color: '#17A2B8', icon: 'fas fa-arrow-up' },
    'Vorticidade Total (1/s)': { label: 'Vorticidade Absoluta Total', color: '#9C27B0', icon: 'fas fa-hurricane' },
    'Pressão Média (Pa)': { label: 'Pressão Média', color: '#fde725', icon: 'fas fa-gauge-high' },
    'Número CFL Máx': { label: 'Número CFL Máximo', color: '#E91E63', icon: 'fas fa-triangle-exclamation'},
    'Velocidade Máx (m/s)': { label: 'Velocidade Máxima', color: '#6610f2', icon: 'fas fa-tachometer-alt'},
    'Residual Div. (L2)': { label: 'Residual da Divergência', color: '#673e00ff', icon: 'fas fa-check-circle' }
};

const STATS_CONFIG = {
    avg: { title: 'Média', icon: `<svg class="icon" viewBox="0 0 24 24" style="width:16px; height:16px;"><path fill="currentColor" d="M16 11.78l4.24-7.33l-1.73-1l-3.36 5.82l-2.43-4.22l-1.73 1l1.58 2.74l-4.18 7.23l-1.73-1l-3.36 5.82l1.73 1l2.48-4.3l2.43 4.22l1.73-1l-1.58-2.74l4.18-7.23z"/></svg>` },
    max: { title: 'Pico (Máx)', icon: `<svg class="icon" viewBox="0 0 24 24" style="width:16px; height:16px;"><path fill="currentColor" d="M16 6l2.29 2.29l-4.88 4.88l-4-4L2 16.59L3.41 18l6-6l4 4l6.3-6.29L22 12V6z"/></svg>` },
    median: { title: 'Mediana', icon: `<svg class="icon" viewBox="0 0 24 24" style="width:16px; height:16px;"><path fill="currentColor" d="M4 18h3V6H4v12M8 18h3V9H8v9m4 0h3v-5h-3v5m4 0h3V3h-3v15z"/></svg>`},
    stdDev: { title: 'Desv. Padrão', icon: `<svg class="icon" viewBox="0 0 24 24" style="width:16px; height:16px;"><path fill="currentColor" d="M3.5 18.5c.83 0 1.5-.67 1.5-1.5V7H7v10c0 2.21 1.79 4 4 4s4-1.79 4-4V7h1.5v10c0 .83.67 1.5 1.5 1.5s1.5-.67 1.5-1.5V6c0-1.11-.89-2-2-2h-1c-1.11 0-2 .89-2 2v1h-1.5V6c0-1.11-.89-2-2-2s-2 .89-2 2v1H9V6c0-1.11-.89-2-2-2H6c-1.11 0-2 .89-2 2v11.5c0 .83.67 1.5 1.5 1.5z"/></svg>`},
    integral: { title: 'Integral (·s)', icon: `<svg class="icon" viewBox="0 0 24 24" style="width:16px; height:16px;"><path fill="currentColor" d="M2 20h20V4l-5.32 7.1l-4.69-3.52L6 14.24l-4 3.73V20m2-11l3-2.73l5.09 3.82l4.69-3.52L20 9.7V6H4v5.36l-2 1.87V6h0z"/></svg>`},
};

function populateReportModal() {
    correlationChartsRendered = false;

    const paramsList = document.getElementById('paramsList');
    paramsList.innerHTML = '';
    for (const [key, value] of Object.entries(simulationParametersSnapshot)) {
        const formattedValue = typeof value === 'number' ? formatNumberForDisplay(value) : value;
        paramsList.innerHTML += `<li>${key}: <b>${formattedValue}</b></li>`;
    }
    const duration = recordedData.length > 0 ? recordedData[recordedData.length - 1]['Tempo (s)'] : 0;
    paramsList.innerHTML += `<li>Duração da Gravação: <b>${duration.toFixed(1)}s</b></li>`;
    paramsList.innerHTML += `<li>Frequência de Amostragem: <b>${simulationParametersSnapshot['Frequência de Amostragem (Hz)'] || recordingFrequency} Hz</b></li>`;
    paramsList.innerHTML += `<li>Total de Amostras: <b>${recordedData.length}</b></li>`;

    const statsList = document.getElementById('statsList');
    statsList.innerHTML = '';
    for(const key of Object.keys(METRIC_CONFIG)) {
        if (!recordedData[0] || !recordedData[0].hasOwnProperty(key)) continue;
        const stats = getAggregateStats(recordedData, key); const config = METRIC_CONFIG[key];
        let subListHtml = `<ul class="sub-stats-list">`;
        for (const statKey of Object.keys(STATS_CONFIG)) {
            const statConfig = STATS_CONFIG[statKey]; const statValue = stats[statKey];
            // --- INÍCIO DA REFATORAÇÃO: Usar nova função de formatação para as estatísticas ---
            subListHtml += `<li>${statConfig.icon} ${statConfig.title}: <b class="monospace-font">${formatNumberForDisplay(statValue, 3)}</b></li>`;
            // --- FIM DA REFATORAÇÃO ---
        }
        subListHtml += `</ul>`;
        statsList.innerHTML += `<li><div class="metric-header"><i class="${config.icon}" style="color: ${config.color}"></i> ${config.label}</div>${subListHtml}</li>`;
    }
    populateTable(recordedData); 
    setupTableSorting(); 
    drawAllCharts(); 
    setupCustomChart();
}

function drawAllCharts() {
    const chartsGrid = document.getElementById('chartsGrid');
    chartsGrid.innerHTML = '';
    const useLogScale = document.getElementById('logScaleToggle').checked;
    const showMovingAverage = document.getElementById('movingAverageToggle').checked;
    for(const key of Object.keys(METRIC_CONFIG)){
        if (!recordedData[0] || !recordedData[0].hasOwnProperty(key)) continue;
        const config = METRIC_CONFIG[key]; const wrapper = document.createElement('div');
        wrapper.className = 'chart-wrapper'; const canvas = document.createElement('canvas');
        wrapper.appendChild(canvas); chartsGrid.appendChild(wrapper);
        const movingAverageData = showMovingAverage ? calculateMovingAverage(recordedData.map(d => d[key]), 5) : null;
        let hoveredIndex = -1;
        const draw = () => { drawChart(canvas, recordedData, key, config.label, config.color, hoveredIndex, useLogScale, movingAverageData); };
        wrapper.addEventListener('mousemove', (e) => {
            const rect = canvas.getBoundingClientRect(); const x = e.clientX - rect.left;
            const padding = { left: 50, right: 20 }; const chartW = rect.width - padding.left - padding.right;
            let newIndex = Math.round(((x) / chartW) * (recordedData.length - 1));
            newIndex = Math.max(0, Math.min(recordedData.length - 1, newIndex));
            if (newIndex !== hoveredIndex) { hoveredIndex = newIndex; draw(); }
        });
        wrapper.addEventListener('mouseleave', () => { hoveredIndex = -1; draw(); });
        draw();
    }
}
document.getElementById('logScaleToggle').addEventListener('change', drawAllCharts);
document.getElementById('movingAverageToggle').addEventListener('change', drawAllCharts);
function populateTable(data) {
    const tableHead = document.querySelector('#reportTable thead');
    const tableBody = document.querySelector('#reportTable tbody');
    tableHead.innerHTML = ''; tableBody.innerHTML = '';
    if (data.length === 0) return;
    const filteredData = data.map(row => { const newRow = { ...row }; delete newRow['Histograma EC']; return newRow; });
    const headers = Object.keys(filteredData[0]);
    let headerHtml = '<tr>';
    headers.forEach(h => headerHtml += `<th class="sortable" data-key="${h}">${h}<span class="sort-arrow"></span></th>`);
    headerHtml += '</tr>'; tableHead.innerHTML = headerHtml;
    let tableHtml = '';
    filteredData.forEach(row => {
        tableHtml += '<tr>';
        // --- INÍCIO DA REFATORAÇÃO: Formatar números na tabela ---
        headers.forEach(h => tableHtml += `<td>${typeof row[h] === 'number' ? formatNumberForDisplay(row[h], 3) : row[h]}</td>`);
        // --- FIM DA REFATORAÇÃO ---
        tableHtml += '</tr>';
    });
    tableBody.innerHTML = tableHtml;
}
let sortState = { key: null, asc: true };
function setupTableSorting() {
    document.querySelectorAll('#reportTable th.sortable').forEach(th => {
        th.addEventListener('click', () => {
            const key = th.dataset.key;
            if (sortState.key === key) { sortState.asc = !sortState.asc; } else { sortState.key = key; sortState.asc = true; }
            recordedData.sort((a, b) => {
                if (a[key] < b[key]) return sortState.asc ? -1 : 1;
                if (a[key] > b[key]) return sortState.asc ? 1 : -1;
                return 0;
            });
            populateTable(recordedData);
            document.querySelectorAll('#reportTable th .sort-arrow').forEach(arrow => arrow.textContent = '');
            th.querySelector('.sort-arrow').textContent = sortState.asc ? '▲' : '▼';
        });
    });
}
tableSearchInput.addEventListener('input', (e) => {
    const searchTerm = e.target.value.toLowerCase();
    if (!searchTerm) { recordedData = [...originalRecordedData]; populateTable(recordedData); return; }
    recordedData = originalRecordedData.filter(row => Object.values(row).some(value => value.toString().toLowerCase().includes(searchTerm)));
    populateTable(recordedData);
});
function calculateMovingAverage(data, windowSize) {
    if (windowSize > data.length) return []; let result = [];
    for (let i = 0; i < data.length; i++) {
        const start = Math.max(0, i - windowSize + 1); const window = data.slice(start, i + 1);
        const avg = window.reduce((a, b) => a + b, 0) / window.length;
        result.push(avg);
    }
    return result;
}
function drawChart(canvas, data, dataKey, label, color, hoveredIndex = -1, useLogScale = false, trendData = null) {
    if (!canvas || !data || data.length === 0 || !data[0].hasOwnProperty(dataKey)) return;
    const dpr = window.devicePixelRatio || 1; const rect = canvas.getBoundingClientRect();
    canvas.width = rect.width * dpr; canvas.height = rect.height * dpr;
    const ctx = canvas.getContext('2d'); ctx.scale(dpr, dpr);
    const w = rect.width; const h = rect.height; ctx.clearRect(0, 0, w, h);
    const values = data.map(d => d[dataKey]);
    const safeValues = useLogScale ? values.map(v => v > 1e-12 ? v : 1e-12) : values;
    let maxVal = Math.max(...safeValues); let minVal = Math.min(...safeValues);
    if (useLogScale && maxVal <= 0) { maxVal = 1; minVal = 1e-12; }
    if (!useLogScale && maxVal > 0 && minVal > maxVal * 0.95) minVal = 0;
    const padding = { top: 30, bottom: 20, left: 50, right: 20 };
    const chartW = w - padding.left - padding.right; const chartH = h - padding.top - padding.bottom;
    const mapY = (val) => {
        if (useLogScale) {
            const logMin = Math.log10(minVal); const logMax = Math.log10(maxVal);
            const logRange = logMax - logMin < 1e-9 ? 1 : logMax - logMin;
            return padding.top + chartH - ((Math.log10(val) - logMin) / logRange) * chartH;
        } else {
            const range = maxVal - minVal < 1e-9 ? 1 : maxVal - minVal;
            return padding.top + chartH - ((val - minVal) / range) * chartH;
        }
    };
    ctx.strokeStyle = '#333'; ctx.lineWidth = 0.5; ctx.font = "10px 'Roboto Mono'"; ctx.fillStyle = 'var(--text-secondary)'; ctx.textAlign = 'right';
    const numGridLines = 4;
    for(let i=0; i <= numGridLines; i++){
        const y = padding.top + (chartH / numGridLines) * i;
        ctx.beginPath(); ctx.moveTo(padding.left, y); ctx.lineTo(padding.left + chartW, y); ctx.stroke();
        let val;
        if(useLogScale) { const logMin = Math.log10(minVal), logMax = Math.log10(maxVal); val = Math.pow(10, logMax - ((logMax - logMin) / numGridLines * i)); }
        else { val = maxVal - ((maxVal - minVal) / numGridLines * i); }
        // Nota: A renderização no canvas não suporta HTML, então a notação 'e' é mantida aqui.
        ctx.fillText(val.toExponential(1), padding.left - 8, y + 3);
    }
    ctx.strokeStyle = color; ctx.lineWidth = 1.5; ctx.beginPath();
    for (let i = 0; i < data.length; i++) {
        const x = padding.left + (i / (data.length - 1)) * chartW; const y = mapY(safeValues[i]);
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    if(trendData) {
        ctx.strokeStyle = 'rgba(104, 104, 104, 0.4)'; ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath();
        for (let i = 0; i < trendData.length; i++) {
            const x = padding.left + (i / (trendData.length - 1)) * chartW; const y = mapY(trendData[i] > 0 ? trendData[i] : (useLogScale ? 1e-12 : 0));
            if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        ctx.stroke(); ctx.setLineDash([]);
    }
    ctx.fillStyle = 'var(--text-primary)'; ctx.font = "12px 'Inter'"; ctx.textAlign = 'left'; ctx.fillText(label, padding.left, padding.top - 10);
    if (hoveredIndex !== -1 && data[hoveredIndex]) {
        const point = data[hoveredIndex]; const pointValue = safeValues[hoveredIndex];
        const x = padding.left + (hoveredIndex / (data.length - 1)) * chartW; const y = mapY(pointValue);
        ctx.strokeStyle = 'rgba(136, 136, 136, 0.5)'; ctx.lineWidth = 0.5;
        ctx.beginPath(); ctx.moveTo(x, padding.top); ctx.lineTo(x, padding.top + chartH); ctx.stroke();
        ctx.beginPath(); ctx.arc(x, y, 4, 0, 2 * Math.PI); ctx.fillStyle = color; ctx.fill(); ctx.strokeStyle = '#000'; ctx.lineWidth = 1.5; ctx.stroke();
        const lines = [`Tempo: ${point['Tempo (s)']}s`, `${dataKey}: ${point[dataKey].toExponential(3)}`];
        ctx.font = "11px 'Roboto Mono'";
        const boxWidth = Math.max(...lines.map(t => ctx.measureText(t).width)) + 16; const boxHeight = lines.length * 14 + 8;
        let boxX = x + 15, boxY = y - boxHeight / 2;
        if (boxX + boxWidth > w) boxX = x - 15 - boxWidth;
        boxY = Math.max(padding.top, Math.min(boxY, h - padding.bottom - boxHeight));
        ctx.fillStyle = 'rgba(30, 30, 30, 0.9)'; ctx.strokeStyle = color; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.roundRect(boxX, boxY, boxWidth, boxHeight, 4); ctx.fill(); ctx.stroke();
        ctx.fillStyle = '#ffffff'; ctx.textAlign = 'left';
        lines.forEach((t, i) => ctx.fillText(t, boxX + 8, boxY + 16 + i * 14));
    }
}
function setupCustomChart() {
    const controlsContainer = document.getElementById('custom-chart-controls');
    const canvas = document.getElementById('customChartCanvas'); let hoveredIndex = -1;
    controlsContainer.innerHTML = '';
    Object.keys(METRIC_CONFIG).forEach((metricKey, index) => {
        if (!recordedData[0] || !recordedData[0].hasOwnProperty(metricKey)) return;
        const config = METRIC_CONFIG[metricKey]; const isChecked = ['Energia Cinética (J)', 'Energia Térmica (J)'].includes(metricKey);
        controlsContainer.innerHTML += `<div class="checkbox-container"><input id="custom-check-${metricKey}" type="checkbox" data-metric="${metricKey}" ${isChecked ? 'checked' : ''}/><label for="custom-check-${metricKey}" style="color: ${config.color}">${metricKey}</label></div>`;
    });
    controlsContainer.innerHTML += `<div class="checkbox-container" style="margin-left: auto;"><input id="custom-normalize" type="checkbox"/><label for="custom-normalize">Normalizar Eixos (0-1)</label></div>`;
    controlsContainer.querySelectorAll('input[type="checkbox"]').forEach(checkbox => {
        checkbox.addEventListener('change', () => drawCustomChart(canvas, hoveredIndex));
    });
    const wrapper = document.getElementById('custom-chart-wrapper');
    wrapper.addEventListener('mousemove', (e) => {
        const rect = canvas.getBoundingClientRect(); const x = e.clientX - rect.left;
        const padding = { left: 50, right: 20 }; const chartW = rect.width - padding.left - padding.right;
        let newIndex = Math.round((x / chartW) * (recordedData.length - 1));
        newIndex = Math.max(0, Math.min(recordedData.length - 1, newIndex));
        if (newIndex !== hoveredIndex) { hoveredIndex = newIndex; drawCustomChart(canvas, hoveredIndex); }
    });
    wrapper.addEventListener('mouseleave', () => { hoveredIndex = -1; drawCustomChart(canvas, -1); });
    drawCustomChart(canvas, -1);
}
function drawCustomChart(canvas, hoveredIndex) {
    const controls = document.getElementById('custom-chart-controls');
    const normalize = controls.querySelector('#custom-normalize').checked;
    const selectedMetrics = Array.from(controls.querySelectorAll('input[data-metric]:checked')).map(cb => cb.dataset.metric);
    const dpr = window.devicePixelRatio || 1; const rect = canvas.getBoundingClientRect();
    canvas.width = rect.width * dpr; canvas.height = rect.height * dpr;
    const ctx = canvas.getContext('2d'); ctx.scale(dpr, dpr);
    const w = rect.width; const h = rect.height; ctx.clearRect(0, 0, w, h);
    if (selectedMetrics.length === 0) { ctx.fillStyle = 'var(--text-secondary)'; ctx.font = "14px 'Inter'"; ctx.textAlign = 'center'; ctx.fillText('Selecione uma ou mais métricas para comparar.', w / 2, h / 2); return; }
    const datasets = {}; let globalMin = Infinity, globalMax = -Infinity;
    selectedMetrics.forEach(key => {
        if (recordedData.length > 0 && recordedData[0].hasOwnProperty(key)) {
            const values = recordedData.map(d => d[key]); const min = Math.min(...values); const max = Math.max(...values);
            const range = max - min < 1e-9 ? 1 : max - min;
            if (normalize) { datasets[key] = values.map(v => (v - min) / range); } else { datasets[key] = values; globalMin = Math.min(globalMin, min); globalMax = Math.max(globalMax, max); }
        }
    });
    const yMin = normalize ? 0 : globalMin; const yMax = normalize ? 1 : globalMax;
    const yRange = yMax - yMin < 1e-9 ? 1 : yMax - yMin;
    const padding = { top: 30, bottom: 20, left: 50, right: 20 };
    const chartW = w - padding.left - padding.right; const chartH = h - padding.top - padding.bottom;
    const mapY = (val) => padding.top + chartH - ((val - yMin) / yRange) * chartH;
    ctx.strokeStyle = '#333'; ctx.lineWidth = 0.5; ctx.font = "10px 'Roboto Mono'"; ctx.fillStyle = 'var(--text-secondary)'; ctx.textAlign = 'right';
    const numGridLines = 4;
    for (let i = 0; i <= numGridLines; i++) {
        const y = padding.top + (chartH / numGridLines) * i;
        ctx.beginPath(); ctx.moveTo(padding.left, y); ctx.lineTo(padding.left + chartW, y); ctx.stroke();
        const val = yMax - (yRange / numGridLines * i);
        ctx.fillText(normalize ? val.toFixed(2) : val.toExponential(1), padding.left - 8, y + 3);
    }
    selectedMetrics.forEach(key => {
        if(datasets[key]) {
            const config = METRIC_CONFIG[key]; const values = datasets[key];
            ctx.strokeStyle = config.color; ctx.lineWidth = 1.5; ctx.beginPath();
            for (let i = 0; i < values.length; i++) {
                const x = padding.left + (i / (values.length - 1)) * chartW; const y = mapY(values[i]);
                if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
            }
            ctx.stroke();
        }
    });
    ctx.font = "11px 'Inter'"; let legendX = padding.left;
    selectedMetrics.forEach(key => {
        const config = METRIC_CONFIG[key]; ctx.fillStyle = config.color;
        ctx.fillText('■', legendX, padding.top - 10); ctx.fillStyle = 'var(--text-primary)';
        ctx.textAlign = 'left'; const textWidth = ctx.measureText(config.label).width;
        ctx.fillText(config.label, legendX + 12, padding.top - 10);
        legendX += textWidth + 30;
    });
    if (hoveredIndex !== -1 && recordedData[hoveredIndex]) {
        const point = recordedData[hoveredIndex]; const x = padding.left + (hoveredIndex / (recordedData.length - 1)) * chartW;
        ctx.strokeStyle = 'rgba(136, 136, 136, 0.5)'; ctx.lineWidth = 0.5;
        ctx.beginPath(); ctx.moveTo(x, padding.top); ctx.lineTo(x, padding.top + chartH); ctx.stroke();
        selectedMetrics.forEach(key => {
            if(datasets[key]) {
                const val = datasets[key][hoveredIndex]; const y = mapY(val);
                ctx.beginPath(); ctx.arc(x, y, 4, 0, 2 * Math.PI); ctx.fillStyle = METRIC_CONFIG[key].color; ctx.fill();
                ctx.strokeStyle = '#000'; ctx.lineWidth = 1.5; ctx.stroke();
            }
        });
        let lines = [`Tempo: ${point['Tempo (s)']}s`];
        selectedMetrics.forEach(key => {
            if (point.hasOwnProperty(key)) lines.push(`${key}: ${point[key].toExponential(3)}`);
        });
        ctx.font = "11px 'Roboto Mono'";
        const boxWidth = Math.max(...lines.map(t => ctx.measureText(t).width)) + 16; const boxHeight = lines.length * 14 + 8;
        let boxX = x + 15, boxY = h/3;
        if (boxX + boxWidth > w) boxX = x - 15 - boxWidth;
        ctx.fillStyle = 'rgba(30, 30, 30, 0.9)'; ctx.strokeStyle = 'var(--border-color)'; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.roundRect(boxX, boxY, boxWidth, boxHeight, 4); ctx.fill(); ctx.stroke();
        ctx.textAlign = 'left';
        lines.forEach((t, i) => { const color = i > 0 ? METRIC_CONFIG[selectedMetrics[i-1]].color : '#ffffff'; ctx.fillStyle = color; ctx.fillText(t, boxX + 8, boxY + 16 + i * 14); });
    }
}

function populateCsvOptions() {
    const container = document.getElementById('csv-options-body');
    if (container.childElementCount > 0 || recordedData.length === 0) return;

    const headers = Object.keys(recordedData[0]).filter(h => h !== 'Histograma EC');
    headers.forEach(header => {
        const div = document.createElement('div');
        div.className = 'checkbox-container';
        div.innerHTML = `
            <input id="csv-col-${header}" type="checkbox" data-header="${header}" checked>
            <label for="csv-col-${header}">${header}</label>
        `;
        container.appendChild(div);
    });
}

function handleCsvDownload() {
    const selectedHeaders = Array.from(document.querySelectorAll('#csv-options-body input[type=checkbox]:checked'))
                                 .map(cb => cb.dataset.header);

    if (selectedHeaders.length === 0) {
        alert("Por favor, selecione pelo menos uma coluna para baixar.");
        return;
    }
    
    downloadCSV(selectedHeaders);
    csvOptionsModal.style.display = 'none';
}

function downloadCSV(headers) {
    if (recordedData.length === 0) return;
    
    const csvData = recordedData.map(row => {
        const newRow = {};
        headers.forEach(header => {
            newRow[header] = row[header];
        });
        return newRow;
    });
    
    let csvContent = "data:text/csv;charset=utf-8," + headers.join(',') + '\n';
    csvData.forEach(row => {
        csvContent += headers.map(header => row[header]).join(',') + '\n';
    });
    
    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", `sim_report_${new Date().toISOString()}.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

function handlePdfDownload() {
    const sections = {
        summary: document.getElementById('pdf-option-summary').checked,
        charts: document.getElementById('pdf-option-charts').checked,
        data: document.getElementById('pdf-option-data').checked,
    };
    
    const summaryEl = document.getElementById('tab-summary');
    const chartsEl = document.getElementById('tab-charts');
    const dataEl = document.getElementById('tab-data');
    
    if (!sections.summary) summaryEl.classList.add('section-hidden-in-print');
    if (!sections.charts) chartsEl.classList.add('section-hidden-in-print');
    if (!sections.data) dataEl.classList.add('section-hidden-in-print');
    
    window.print();
    pdfOptionsModal.style.display = 'none';

    setTimeout(() => {
        summaryEl.classList.remove('section-hidden-in-print');
        chartsEl.classList.remove('section-hidden-in-print');
        dataEl.classList.remove('section-hidden-in-print');
    }, 500);
}

// --- Event Listeners dos Modais ---
recordBtn.addEventListener('click', startRecording);
closeReportBtn.addEventListener('click', () => reportModal.style.display = 'none');

document.querySelectorAll('.tab-link').forEach(tab => {
    tab.addEventListener('click', e => {
        const targetTab = e.target.dataset.tab;
        
        document.querySelectorAll('.tab-link').forEach(t => t.classList.remove('active'));
        document.querySelectorAll('.modal-tab-content').forEach(content => content.classList.remove('active'));

        e.target.classList.add('active');
        document.getElementById(`tab-${targetTab}`).classList.add('active');

        // --- INÍCIO DA CORREÇÃO ---
        // Se a aba de gráficos for a selecionada, redesenha todos os gráficos
        // para garantir que sejam renderizados corretamente, mesmo que a aba
        // estivesse oculta anteriormente.
        if (targetTab === 'charts') {
            drawAllCharts();
        }
        // --- FIM DA CORREÇÃO ---

        if (targetTab === 'correlations' && !correlationChartsRendered) {
            populateCorrelationsTab();
            correlationChartsRendered = true;
        }
    });
});

downloadPdfBtn.addEventListener('click', () => pdfOptionsModal.style.display = 'flex');
downloadCsvBtn.addEventListener('click', () => {
    populateCsvOptions();
    csvOptionsModal.style.display = 'flex';
});

document.getElementById('confirmPdfDownloadBtn').addEventListener('click', handlePdfDownload);
document.getElementById('closePdfOptionsBtn').addEventListener('click', () => pdfOptionsModal.style.display = 'none');
document.getElementById('cancelPdfOptionsBtn').addEventListener('click', () => pdfOptionsModal.style.display = 'none');

document.getElementById('confirmCsvDownloadBtn').addEventListener('click', handleCsvDownload);
document.getElementById('closeCsvOptionsBtn').addEventListener('click', () => csvOptionsModal.style.display = 'none');
document.getElementById('cancelCsvOptionsBtn').addEventListener('click', () => csvOptionsModal.style.display = 'none');


const styleSheet = document.createElement("style");
styleSheet.innerText = `@keyframes pulse { 0% { opacity: 1; } 50% { opacity: 0.3; } 100% { opacity: 1; } }`;
document.head.appendChild(styleSheet);


// ===================================================================
// --- INÍCIO: NOVAS FUNÇÕES PARA A ABA DE CORRELAÇÕES ---
// ===================================================================

// --- FUNÇÕES DE CÁLCULO ESTATÍSTICO ---

function calculateRateOfChange(data, key, timeKey = 'Tempo (s)') {
    const derivatives = [0];
    for (let i = 1; i < data.length; i++) {
        const d_val = data[i][key] - data[i-1][key];
        const d_time = data[i][timeKey] - data[i-1][timeKey];
        derivatives.push(d_time > 1e-9 ? d_val / d_time : 0);
    }
    return derivatives;
}

function calculateCorrelation(data, keyX, keyY) {
    if (data.length < 2) return 0;
    let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
    for (const d of data) {
        const x = d[keyX], y = d[keyY];
        sumX += x; sumY += y; sumXY += x * y; sumX2 += x * x; sumY2 += y * y;
    }
    const n = data.length;
    const numerator = n * sumXY - sumX * sumY;
    const denominator = Math.sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));
    return denominator < 1e-9 ? 0 : numerator / denominator;
}

function calculateLinearRegression(data, keyX, keyY) {
    if (data.length < 2) return { m: 0, b: 0, r2: 0 };
    let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    for (const d of data) {
        const x = d[keyX], y = d[keyY];
        sumX += x; sumY += y; sumXY += x * y; sumX2 += x * x;
    }
    const n = data.length;
    const m = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    const b = (sumY - m * sumX) / n;
    const r = calculateCorrelation(data, keyX, keyY);
    return { m: m || 0, b: b || 0, r2: r * r };
}

function calculateRiseTime(data, key) {
    if (data.length < 2) return { time: 0, start: 0, end: 0 };
    const values = data.map(d => d[key]);
    const maxVal = Math.max(...values);
    const target10 = maxVal * 0.1;
    const target90 = maxVal * 0.9;
    let time10 = -1, time90 = -1;
    for(let i=0; i<data.length; i++) {
        if(time10 === -1 && values[i] >= target10) time10 = data[i]['Tempo (s)'];
        if(time90 === -1 && values[i] >= target90) time90 = data[i]['Tempo (s)'];
    }
    if(time10 !== -1 && time90 !== -1) {
        return { time: time90 - time10, start: time10, end: time90 };
    }
    return { time: 0, start: 0, end: 0 };
}

function calculateOvershoot(data, key) {
    if (data.length < 10) return { peak: 0, steady: 0, overshoot: 0 };
    const values = data.map(d => d[key]);
    const peak = Math.max(...values);
    const last20percent = values.slice(Math.floor(values.length * 0.8));
    const steadyState = last20percent.reduce((a,b) => a+b, 0) / last20percent.length;
    const overshoot = steadyState > 1e-9 ? ((peak - steadyState) / steadyState) * 100 : 0;
    return { peak, steady: steadyState, overshoot };
}

function calculateDistributionStats(histData) {
    const { histogram, binSize, totalCells } = histData;
    if (!totalCells) return { mean: 0, stdDev: 0, skewness: 0, kurtosis: 0 };
    let sum = 0, sum2 = 0, sum3 = 0, sum4 = 0;
    
    histogram.forEach((count, i) => {
        const value = (i + 0.5) * binSize;
        sum += value * count;
        sum2 += value * value * count;
    });

    const mean = sum / totalCells;
    const variance = (sum2 / totalCells) - (mean * mean);
    const stdDev = Math.sqrt(variance);

    if (stdDev < 1e-9) return { mean, stdDev, skewness: 0, kurtosis: 0 };
    
    histogram.forEach((count, i) => {
        const value = (i + 0.5) * binSize;
        sum3 += Math.pow(value - mean, 3) * count;
        sum4 += Math.pow(value - mean, 4) * count;
    });

    const skewness = sum3 / (totalCells * Math.pow(stdDev, 3));
    const kurtosis = (sum4 / (totalCells * Math.pow(stdDev, 4))) - 3; // Excesso de curtose

    return { mean, stdDev, skewness, kurtosis };
}


// --- FUNÇÕES DE RENDERIZAÇÃO DE GRÁFICOS (ATUALIZADAS COM TOOLTIPS) ---
function drawScatterPlot(canvas, data, xKey, yKey, title, color, regression = null, hoveredIndex = -1) {
    if (!canvas || data.length < 2) return;
    const dpr = window.devicePixelRatio || 1; const rect = canvas.getBoundingClientRect();
    canvas.width = rect.width * dpr; canvas.height = rect.height * dpr;
    const ctx = canvas.getContext('2d'); ctx.scale(dpr, dpr);
    const w = rect.width; const h = rect.height; ctx.clearRect(0, 0, w, h);

    // Aumentamos o padding para dar espaço para os novos rótulos e títulos
    const padding = { top: 40, bottom: 50, left: 70, right: 20 };

    const xVals = data.map(d => d[xKey]); const yVals = data.map(d => d[yKey]);
    let xMin = Math.min(...xVals), xMax = Math.max(...xVals);
    let yMin = Math.min(...yVals), yMax = Math.max(...yVals);
    
    const xRange = xMax - xMin < 1e-9 ? 1 : xMax - xMin;
    const yRange = yMax - yMin < 1e-9 ? 1 : yMax - yMin;

    // Adiciona uma pequena margem para os pontos não ficarem colados nas bordas
    xMin -= xRange * 0.05; xMax += xRange * 0.05; 
    yMin -= yRange * 0.05; yMax += yRange * 0.05;
    
    const chartW = w - padding.left - padding.right; 
    const chartH = h - padding.top - padding.bottom;
    const mapX = val => padding.left + ((val - xMin) / (xMax - xMin)) * chartW;
    const mapY = val => padding.top + chartH * (1 - (val - yMin) / (yMax - yMin));

    // --- INÍCIO DA ADIÇÃO: GRADE E RÓTULOS ---

    const numGridLines = 4;
    ctx.strokeStyle = '#333'; // Cor sutil para a grade
    ctx.lineWidth = 0.5;
    ctx.font = "10px 'Roboto Mono'";
    ctx.fillStyle = 'var(--text-secondary)';

    // Desenha as linhas da grade horizontal e seus rótulos no eixo Y
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let i = 0; i <= numGridLines; i++) {
        const value = yMin + (yRange / numGridLines) * i;
        const y = mapY(value);
        ctx.beginPath();
        ctx.moveTo(padding.left, y);
        ctx.lineTo(padding.left + chartW, y);
        ctx.stroke();
        ctx.fillText(value.toExponential(1), padding.left - 8, y);
    }

    // Desenha as linhas da grade vertical e seus rótulos no eixo X
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (let i = 0; i <= numGridLines; i++) {
        const value = xMin + (xRange / numGridLines) * i;
        const x = mapX(value);
        ctx.beginPath();
        ctx.moveTo(x, padding.top);
        ctx.lineTo(x, padding.top + chartH);
        ctx.stroke();
        ctx.fillText(value.toExponential(1), x, padding.top + chartH + 8);
    }
    
    // Adiciona os TÍTULOS dos eixos
    ctx.fillStyle = 'var(--text-primary)';
    ctx.font = "12px 'Inter'";
    // Título do Eixo X
    ctx.fillText(xKey, padding.left + chartW / 2, padding.top + chartH + 25);
    
    // Título do Eixo Y (desenhado verticalmente)
    ctx.save();
    ctx.translate(padding.left - 50, padding.top + chartH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(yKey, 0, 0);
    ctx.restore();

    // --- FIM DA ADIÇÃO ---

    // Desenha os pontos de dados
    ctx.fillStyle = color;
    for(let i=0; i<data.length; i++) {
        ctx.beginPath();
        ctx.arc(mapX(xVals[i]), mapY(yVals[i]), 2, 0, 2 * Math.PI);
        ctx.fill();
    }
    
    // Desenha a linha de regressão (se houver)
    if (regression) {
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.7)'; ctx.lineWidth = 1.5; ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(mapX(xMin), mapY(regression.m * xMin + regression.b));
        ctx.lineTo(mapX(xMax), mapY(regression.m * xMax + regression.b));
        ctx.stroke();
        ctx.setLineDash([]);
    }
    
    // Título principal do gráfico
    ctx.fillStyle = 'var(--text-primary)'; ctx.font = "bold 14px 'Inter'"; ctx.textAlign = 'left';
    ctx.fillText(title, padding.left, padding.top - 15);
    
    // Lógica do tooltip (inalterada)
    if (hoveredIndex !== -1 && data[hoveredIndex]) {
        const point = data[hoveredIndex];
        const xVal = point[xKey], yVal = point[yKey];
        const x = mapX(xVal), y = mapY(yVal);

        ctx.strokeStyle = 'rgba(136, 136, 136, 0.5)'; ctx.lineWidth = 0.5;
        ctx.beginPath(); ctx.moveTo(x, padding.top); ctx.lineTo(x, padding.top + chartH); ctx.stroke();
        
        ctx.beginPath(); ctx.arc(x, y, 4, 0, 2 * Math.PI); ctx.fillStyle = color; ctx.fill();
        ctx.strokeStyle = '#000'; ctx.lineWidth = 1.5; ctx.stroke();
        
        const lines = [`${xKey.split(' ')[0]}: ${xVal.toExponential(3)}`, `${yKey.split(' ')[0]}: ${yVal.toExponential(3)}`];
        ctx.font = "11px 'Roboto Mono'";
        const boxWidth = Math.max(...lines.map(t => ctx.measureText(t).width)) + 16;
        const boxHeight = lines.length * 14 + 8;
        let boxX = x + 15, boxY = y - boxHeight / 2;
        if (boxX + boxWidth > w) boxX = x - 15 - boxWidth;
        boxY = Math.max(padding.top, Math.min(boxY, h - padding.bottom - boxHeight));

        ctx.fillStyle = 'rgba(30, 30, 30, 0.9)'; ctx.strokeStyle = color; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.roundRect(boxX, boxY, boxWidth, boxHeight, 4); ctx.fill(); ctx.stroke();
        
        ctx.fillStyle = '#ffffff'; ctx.textAlign = 'left';
        lines.forEach((t, i) => ctx.fillText(t, boxX + 8, boxY + 16 + i * 14));
    }
}

function drawHistogram(canvas, histData, title, color, hoveredIndex = -1) {
    if (!canvas || !histData || !histData.histogram) return;
    const { histogram, binSize } = histData;
    const dpr = window.devicePixelRatio || 1; const rect = canvas.getBoundingClientRect();
    canvas.width = rect.width * dpr; canvas.height = rect.height * dpr;
    const ctx = canvas.getContext('2d'); ctx.scale(dpr, dpr);
    const w = rect.width; const h = rect.height; ctx.clearRect(0, 0, w, h);
    
    const maxCount = Math.max(...histogram);
    const padding = { top: 30, bottom: 20, left: 50, right: 20 };
    const chartW = w - padding.left - padding.right;
    const chartH = h - padding.top - padding.bottom;
    const barWidth = chartW / histogram.length;

    for (let i = 0; i < histogram.length; i++) {
        const barHeight = (histogram[i] / maxCount) * chartH;
        ctx.fillStyle = i === hoveredIndex ? '#fff' : color;
        ctx.fillRect(padding.left + i * barWidth, padding.top + chartH - barHeight, barWidth - 1, barHeight);
    }
    
    ctx.fillStyle = 'var(--text-primary)'; ctx.font = "12px 'Inter'"; ctx.textAlign = 'left';
    ctx.fillText(title, padding.left, padding.top - 10);
    ctx.font = "10px 'Roboto Mono'"; ctx.fillStyle = 'var(--text-secondary)'; ctx.textAlign = 'center';
    ctx.fillText((0).toExponential(1), padding.left, padding.top + chartH + 12);
    ctx.textAlign = 'right';
    ctx.fillText((histogram.length * binSize).toExponential(1), padding.left + chartW, padding.top + chartH + 12);
    
    if(hoveredIndex !== -1 && histogram[hoveredIndex] !== undefined) {
        const count = histogram[hoveredIndex];
        const rangeStart = hoveredIndex * binSize;
        const rangeEnd = (hoveredIndex + 1) * binSize;
        const barX = padding.left + (hoveredIndex + 0.5) * barWidth;
        
        const lines = [`Faixa: ${rangeStart.toExponential(1)} - ${rangeEnd.toExponential(1)}`, `Contagem: ${count}`];
        ctx.font = "11px 'Roboto Mono'";
        const boxWidth = Math.max(...lines.map(t => ctx.measureText(t).width)) + 16;
        const boxHeight = lines.length * 14 + 8;
        
        // --- INÍCIO DA CORREÇÃO ---

        // 1. Calcula a altura e a posição Y do topo da barra específica
        const barHeight = (histogram[hoveredIndex] / maxCount) * chartH;
        const barTopY = padding.top + chartH - barHeight;

        // 2. Calcula a posição X (horizontal) centrada na barra
        let boxX = barX - boxWidth / 2;
        
        // 3. Calcula a posição Y (vertical) para ficar ACIMA da barra
        let boxY = barTopY - boxHeight - 5; // 5px de margem

        // 4. Garante que o tooltip não saia para os lados...
        boxX = Math.max(padding.left, Math.min(boxX, w - padding.right - boxWidth));
        
        // 5. ...e nem para cima!
        if (boxY < padding.top) {
            boxY = barTopY + 5; // Se for sair para cima, posiciona ABAIXO do topo da barra
        }

        // --- FIM DA CORREÇÃO ---
        
        ctx.fillStyle = 'rgba(30, 30, 30, 0.9)'; ctx.strokeStyle = color; ctx.lineWidth = 1;
        ctx.beginPath(); ctx.roundRect(boxX, boxY, boxWidth, boxHeight, 4); ctx.fill(); ctx.stroke();
        
        ctx.fillStyle = '#ffffff'; ctx.textAlign = 'left';
        lines.forEach((t, i) => ctx.fillText(t, boxX + 8, boxY + 16 + i * 14));
    }
}


// --- FUNÇÃO PRINCIPAL DA ABA DE CORRELAÇÕES ---
function populateCorrelationsTab() {
    const container = document.getElementById('analysisGrid');
    if (recordedData.length < 5) {
        container.innerHTML = '<p style="padding: 20px;">Dados insuficientes para análise de correlação. Grave por mais tempo (pelo menos 5 amostras).</p>';
        return;
    }

    let cardsHtml = '';
    
    // Análise 1
    const convectionData = { title: "Eficiência da Convecção", desc: "Mede como a energia térmica (X) é convertida em energia cinética (Y).", xKey: 'Energia Térmica (J)', yKey: 'Energia Cinética (J)', color: '#FF6B35' };
    const corr_conv = calculateCorrelation(recordedData, convectionData.xKey, convectionData.yKey);
    const r2_conv = corr_conv * corr_conv;
    const convectionStatsHtml = `<ul><li>Coeficiente de Correlação (r): <b class="monospace-font">${corr_conv.toFixed(4)}</b></li><li>Coeficiente de Determinação (R²): <b class="monospace-font">${r2_conv.toFixed(4)}</b></li></ul><small>${(r2_conv*100).toFixed(1)}% da variação na energia cinética pode ser explicada pela variação na energia térmica.</small>`;
    cardsHtml += createAnalysisCard(1, convectionData.title, convectionData.desc, convectionStatsHtml);

    // Análise 2
    const dKEdt = calculateRateOfChange(recordedData, 'Energia Cinética (J)');
    const dissipationData = recordedData.map((d, i) => ({ ...d, 'd(KE)/dt': dKEdt[i] }));
    const turbulenceData = { title: "Turbulência e Dissipação", desc: "Relação entre a dissipação de energia (X) e a vorticidade total (Y).", xKey: 'd(KE)/dt', yKey: 'Vorticidade Total (1/s)', color: '#9C27B0' };
    const corr_turb = calculateCorrelation(dissipationData, turbulenceData.xKey, turbulenceData.yKey);
    const sortedByVort = [...dissipationData].sort((a,b) => a[turbulenceData.yKey] - b[turbulenceData.yKey]);
    const q1_len = Math.floor(sortedByVort.length * 0.25);
    const low_vort_slice = sortedByVort.slice(0, q1_len);
    const high_vort_slice = sortedByVort.slice(sortedByVort.length - q1_len);
    const avg_diss_low = low_vort_slice.reduce((s,d) => s + d[turbulenceData.xKey], 0) / low_vort_slice.length;
    const avg_diss_high = high_vort_slice.reduce((s,d) => s + d[turbulenceData.xKey], 0) / high_vort_slice.length;
    const turbulenceStatsHtml = `<ul><li>Correlação (r): <b class="monospace-font">${corr_turb.toFixed(4)}</b></li><li>Dissipação Média (25% menor vort.): <b class="monospace-font">${formatNumberForDisplay(avg_diss_low, 2)} J/s</b></li><li>Dissipação Média (25% maior vort.): <b class="monospace-font">${formatNumberForDisplay(avg_diss_high, 2)} J/s</b></li></ul>`;
    cardsHtml += createAnalysisCard(2, turbulenceData.title, turbulenceData.desc, turbulenceStatsHtml);
    
    // Análise 3
    const stabilityData = { title: "Estabilidade Numérica", desc: "Valida a relação linear entre Velocidade Máxima (X) e Número CFL (Y).", xKey: 'Velocidade Máx (m/s)', yKey: 'Número CFL Máx', color: '#E91E63' };
    const regression = calculateLinearRegression(recordedData, stabilityData.xKey, stabilityData.yKey);
    const stabilityStatsHtml = `<ul><li>Coeficiente de Determinação (R²): <b class="monospace-font">${regression.r2.toFixed(5)}</b></li><li>Equação da Linha: <b class="monospace-font">y = ${formatNumberForDisplay(regression.m, 2)}x + ${formatNumberForDisplay(regression.b, 2)}</b></li></ul><small>Um R² próximo de 1.0 confirma a forte linearidade esperada no solver.</small>`;
    cardsHtml += createAnalysisCard(3, stabilityData.title, stabilityData.desc, stabilityStatsHtml);

    // Análise 4
    const responseData = { title: "Resposta do Sistema à Injeção", desc: "Mede a inércia do sistema, ou seja, quão rápido o fluido reage à fonte de energia.", key: 'Energia Cinética (J)', color: '#007AFF' };
    const riseTime = calculateRiseTime(recordedData, responseData.key);
    const overshoot = calculateOvershoot(recordedData, responseData.key);
    const responseStatsHtml = `<ul><li>Tempo de Subida (10%-90%): <b class="monospace-font">${riseTime.time.toFixed(3)} s</b></li><li>Pico (Overshoot): <b class="monospace-font">${overshoot.overshoot.toFixed(2)} %</b></li><li>Valor de Pico: <b class="monospace-font">${formatNumberForDisplay(overshoot.peak, 2)} J</b></li></ul>`;
    cardsHtml += createAnalysisCard(4, responseData.title, responseData.desc, responseStatsHtml);

    // Análise 5
    const lastFrameHist = recordedData[recordedData.length-1]['Histograma EC'];
    const distData = { title: "Distribuição da Energia Cinética", desc: "Mostra como a energia está distribuída entre as células do fluido no último frame.", color: '#4CAF50'};
    const distStats = calculateDistributionStats(lastFrameHist);
    const distributionStatsHtml = `<ul><li>Assimetria (Skewness): <b class="monospace-font">${distStats.skewness.toFixed(3)}</b></li><li>Curtose (Kurtosis): <b class="monospace-font">${distStats.kurtosis.toFixed(3)}</b></li></ul><small>Assimetria positiva alta indica que a maioria das células tem baixa energia, com uma 'cauda' de alta energia.</small>`;
    cardsHtml += createAnalysisCard(5, distData.title, distData.desc, distributionStatsHtml);
    
    container.innerHTML = cardsHtml;
    
    // --- INÍCIO DA CORREÇÃO ---
    // Função auxiliar para encontrar o ponto mais próximo em um gráfico de dispersão
    const findClosestPointInScatter = (e, canvas, data, xKey, yKey) => {
        const rect = canvas.getBoundingClientRect();
        const mouseX = e.clientX - rect.left;
        const mouseY = e.clientY - rect.top;

        const xVals = data.map(d => d[xKey]);
        const yVals = data.map(d => d[yKey]);
        let xMin = Math.min(...xVals), xMax = Math.max(...xVals);
        let yMin = Math.min(...yVals), yMax = Math.max(...yVals);
        const xRange = xMax - xMin < 1e-9 ? 1 : xMax - xMin;
        const yRange = yMax - yMin < 1e-9 ? 1 : yMax - yMin;
        xMin -= xRange * 0.05; xMax += xRange * 0.05; yMin -= yRange * 0.05; yMax += yRange * 0.05;

        const padding = { top: 30, bottom: 20, left: 50, right: 20 };
        const chartW = rect.width - padding.left - padding.right;
        const chartH = rect.height - padding.top - padding.bottom;
        const mapX = val => padding.left + ((val - xMin) / (xMax - xMin)) * chartW;
        const mapY = val => padding.top + chartH * (1 - (val - yMin) / (yMax - yMin));

        let closestIndex = -1;
        let minDistanceSq = Infinity;

        for (let i = 0; i < data.length; i++) {
            const pointX = mapX(data[i][xKey]);
            const pointY = mapY(data[i][yKey]);
            const distanceSq = (mouseX - pointX)**2 + (mouseY - pointY)**2;
            if (distanceSq < minDistanceSq) {
                minDistanceSq = distanceSq;
                closestIndex = i;
            }
        }
        // Aumenta a área de detecção para ser mais amigável
        return minDistanceSq < 400 ? closestIndex : -1; 
    };
    
    // Gráfico 1: Convecção (Scatter)
    const canvas1 = document.getElementById('corr-chart-1-canvas');
    const wrapper1 = canvas1.parentElement;
    let hoveredIndex1 = -1;
    const draw1 = () => drawScatterPlot(canvas1, recordedData, convectionData.xKey, convectionData.yKey, "Térmica vs. Cinética", convectionData.color, null, hoveredIndex1);
    wrapper1.addEventListener('mousemove', e => {
        const newIndex = findClosestPointInScatter(e, canvas1, recordedData, convectionData.xKey, convectionData.yKey);
        if (newIndex !== hoveredIndex1) {
            hoveredIndex1 = newIndex;
            draw1();
        }
    });
    wrapper1.addEventListener('mouseleave', () => { if (hoveredIndex1 !== -1) { hoveredIndex1 = -1; draw1(); } });
    draw1();

    // Gráfico 2: Turbulência (Scatter)
    const canvas2 = document.getElementById('corr-chart-2-canvas');
    const wrapper2 = canvas2.parentElement;
    let hoveredIndex2 = -1;
    const draw2 = () => drawScatterPlot(canvas2, dissipationData, turbulenceData.xKey, turbulenceData.yKey, "Dissipação vs. Vorticidade", turbulenceData.color, null, hoveredIndex2);
    wrapper2.addEventListener('mousemove', e => {
        const newIndex = findClosestPointInScatter(e, canvas2, dissipationData, turbulenceData.xKey, turbulenceData.yKey);
        if (newIndex !== hoveredIndex2) {
            hoveredIndex2 = newIndex;
            draw2();
        }
    });
    wrapper2.addEventListener('mouseleave', () => { if (hoveredIndex2 !== -1) { hoveredIndex2 = -1; draw2(); } });
    draw2();

    // Gráfico 3: Estabilidade (Scatter)
    const canvas3 = document.getElementById('corr-chart-3-canvas');
    const wrapper3 = canvas3.parentElement;
    let hoveredIndex3 = -1;
    const draw3 = () => drawScatterPlot(canvas3, recordedData, stabilityData.xKey, stabilityData.yKey, "Velocidade Máx vs. CFL", stabilityData.color, regression, hoveredIndex3);
    wrapper3.addEventListener('mousemove', e => {
        const newIndex = findClosestPointInScatter(e, canvas3, recordedData, stabilityData.xKey, stabilityData.yKey);
        if (newIndex !== hoveredIndex3) {
            hoveredIndex3 = newIndex;
            draw3();
        }
    });
    wrapper3.addEventListener('mouseleave', () => { if (hoveredIndex3 !== -1) { hoveredIndex3 = -1; draw3(); } });
    draw3();
    // --- FIM DA CORREÇÃO ---


    // Gráfico 4: Resposta (Time Series) - A lógica original está correta para este tipo de gráfico
    const canvas4 = document.getElementById('corr-chart-4-canvas');
    const wrapper4 = canvas4.parentElement;
    let hoveredIndex4 = -1;
    const draw4 = () => drawChart(canvas4, recordedData, responseData.key, "Energia Cinética vs. Tempo", responseData.color, hoveredIndex4);
    wrapper4.addEventListener('mousemove', (e) => {
        const rect = canvas4.getBoundingClientRect(); const x = e.clientX - rect.left;
        const padding = { left: 50, right: 20 }; const chartW = rect.width - padding.left - padding.right;
        let newIndex = Math.round(((x - padding.left) / chartW) * (recordedData.length - 1));
        newIndex = Math.max(0, Math.min(recordedData.length - 1, newIndex));
        if (newIndex !== hoveredIndex4) { hoveredIndex4 = newIndex; draw4(); }
    });
    wrapper4.addEventListener('mouseleave', () => { if (hoveredIndex4 !== -1) { hoveredIndex4 = -1; draw4(); } });
    draw4();

    // Gráfico 5: Distribuição (Histogram) - A lógica original está correta para este tipo de gráfico
    const canvas5 = document.getElementById('corr-chart-5-canvas');
    const wrapper5 = canvas5.parentElement;
    let hoveredIndex5 = -1;
    const draw5 = () => drawHistogram(canvas5, lastFrameHist, "Histograma de Energia Cinética por Célula", distData.color, hoveredIndex5);
    wrapper5.addEventListener('mousemove', (e) => {
        const rect = canvas5.getBoundingClientRect(); const x = e.clientX - rect.left;
        const padding = { left: 50, right: 20 }; const chartW = rect.width - padding.left - padding.right;
        const barWidth = chartW / lastFrameHist.histogram.length;
        let newIndex = Math.floor((x - padding.left) / barWidth);
        newIndex = Math.max(0, Math.min(lastFrameHist.histogram.length - 1, newIndex));
        if (newIndex !== hoveredIndex5) { hoveredIndex5 = newIndex; draw5(); }
    });
    wrapper5.addEventListener('mouseleave', () => { if (hoveredIndex5 !== -1) { hoveredIndex5 = -1; draw5(); } });
    draw5();
}


function createAnalysisCard(id, title, description, statsHtml) {
    return `
        <div class="analysis-card">
            <h4>${id}. ${title}</h4>
            <p class="description">${description}</p>
            <div class="analysis-content">
                <div class="chart-container">
                    <canvas id="corr-chart-${id}-canvas"></canvas>
                </div>
                <div class="stats-container">${statsHtml}</div>
            </div>
        </div>
    `;
}