# Métodos Numéricos en Crisis

> **Desafío Final — Métodos Numéricos**
> *Simulación numérica de abastecimiento, precios y conflicto social en contexto de crisis*

---

## Descripción

Plataforma web interactiva que aplica métodos numéricos para modelar y simular escenarios reales de crisis: distribución de carburantes, vaciado de reservas, inflación de precios, pérdida del poder adquisitivo familiar y cálculo de umbrales críticos.

Cada escenario permite ingresar parámetros personalizados, ejecutar los cálculos en el navegador y visualizar los resultados mediante tablas, gráficos y explicaciones interpretativas.

---

## Estructura del proyecto

```
proyecto-metodos-numericos/
├── index.html                  ← Dashboard principal
├── css/
│   └── styles.css              ← Estilos compartidos
├── js/
│   └── utils.js                ← Funciones numéricas comunes
└── escenarios/
    ├── escenario-a.html        ← Sistemas de ecuaciones lineales
    ├── escenario-b.html        ← Ecuaciones diferenciales ordinarias
    ├── escenario-c.html        ← Interpolación
    ├── escenario-d.html        ← Integración numérica
    └── escenario-e.html        ← Raíces de ecuaciones
```

---

## Escenarios y métodos implementados

### Escenario A · Optimización de Red de Transporte

Modela la distribución de carburantes y alimentos desde plantas hacia zonas de la ciudad. Permite simular el bloqueo de rutas modificando la matriz de coeficientes.

**Métodos:** Eliminación Gaussiana / LU · Jacobi · Gauss-Seidel · SOR

---

### Escenario B · Vaciado Crítico de Reservas

Simula la evolución de la reserva de carburante con la EDO `R'(t) = entrada − consumo`, con consumo creciente por efecto de pánico.

**Métodos:** Euler · Heun · Runge-Kutta de 4to orden (RK4)

---

### Escenario C · Curva Continua de Precios

Reconstruye la curva de precios de productos básicos a partir de datos dispersos y estima el precio en cualquier día del mes.

**Métodos:** Interpolación de Lagrange · Diferencias divididas de Newton · Splines cúbicos

---

### Escenario D · Gasto Acumulado Familiar

Calcula el costo total que enfrenta una familia durante el mes cuando los precios suben de forma lineal, exponencial o en escalones. El gasto acumulado es el área bajo la curva de precios.

**Métodos:** Regla del Trapecio · Simpson 1/3 · Simpson 3/8

---

### Escenario E · Umbrales Críticos de Abastecimiento

Encuentra el punto exacto donde el sistema "cambia de estado": el día en que el gasto supera el ingreso familiar, la tasa mínima de reposición de carburante, o el umbral de descontento social.

**Métodos:** Bisección · Newton-Raphson · Secante

---

## Cómo ejecutar localmente

1. Clona el repositorio:

   ```bash
   git clone https://github.com/eriekey1/PROYECTO-METODOS-NUMERICOS.git
   ```

2. Abre `index.html` en tu navegador, o usa Live Server en VS Code.

> No requiere instalación de dependencias ni servidor backend. Todo corre en el navegador.

---

## Contexto académico

Este proyecto corresponde al **Desafío Final de Métodos Numéricos**, cuyo objetivo es aplicar los métodos estudiados durante la materia para construir herramientas interactivas capaces de modelar problemas reales. El enfoque es analítico y académico: usar la matemática para entender la realidad, contrastar escenarios y demostrar cómo el modelado numérico apoya la toma de decisiones.

---
## Conclusiones

Este proyecto demostró que los métodos numéricos no son únicamente herramientas abstractas de cálculo, sino instrumentos capaces de modelar fenómenos sociales y económicos reales con precisión y utilidad práctica.

A lo largo de los cinco escenarios desarrollados se pudo comprobar que la eliminación gaussiana y los métodos iterativos permiten representar redes de distribución complejas y simular el efecto de bloqueos sobre el abastecimiento. Los métodos de integración numérica —trapecio, Simpson 1/3 y 3/8— resultaron efectivos para estimar el gasto acumulado familiar bajo distintos patrones de inflación, evidenciando cómo pequeñas variaciones en la curva de precios impactan significativamente el poder adquisitivo. La interpolación de Lagrange, las diferencias divididas de Newton y los splines cúbicos permitieron reconstruir comportamientos continuos a partir de datos dispersos, lo cual es especialmente relevante en contextos donde la información es incompleta o irregular.

La simulación del vaciado de reservas mediante ecuaciones diferenciales ordinarias mostró que el método de Runge-Kutta de cuarto orden ofrece una aproximación notablemente más precisa que Euler o Heun cuando el consumo varía de forma no lineal, como ocurre en situaciones de pánico social. Finalmente, los métodos de búsqueda de raíces —bisección, Newton-Raphson y secante— resultaron ser los más directamente interpretables en términos de política pública, ya que permiten identificar con exactitud el umbral a partir del cual un sistema entra en crisis.
En conjunto, el proyecto confirma que el modelado numérico es una herramienta valiosa para anticipar puntos de quiebre, comparar escenarios y fundamentar decisiones en contextos de incertidumbre. La implementación en el navegador, sin dependencias externas, refuerza además la accesibilidad y replicabilidad de los modelos desarrollados.

---

## Autor

| Campo | Detalle |
|---|---|
| Estudiante | Erika Quecaño Uruña |
| Materia | Métodos Numéricos |
| Universidad | Universidad Mayor de San Andrés |
| Gestión | 2026 |
