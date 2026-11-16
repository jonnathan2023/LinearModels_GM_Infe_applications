# Inferencia en Modelos de Gauss-Markov

[![R](https://img.shields.io/badge/R-%E2%89%A54.0-blue.svg)](https://www.r-project.org/)
[![LaTeX](https://img.shields.io/badge/LaTeX-XeLaTeX-green.svg)](https://www.latex-project.org/)
[![Status](https://img.shields.io/badge/Status-Completado-brightgreen.svg)]()

Repositorio del trabajo encomendado para el curso de **Modelos Lineales I** (Departamento de Estadística e Informática, Universidad Nacional Agraria La Molina - UNALM).

El repositorio aborda la teoría y aplicación de modelos lineales de Gauss-Markov con matrices de diseño de rango deficiente, enfocándose en la estimabilidad de parámetros, distribuciones de formas cuadráticas, intervalos de confianza y pruebas de hipótesis lineales generales.

---

## Contenido Teórico y Práctico

El trabajo resuelve y sistematiza cinco situaciones de inferencia estadística:

1. **Modelación Matricial y Estimabilidad:** Formulación del modelo no de rango completo ($Y = X\beta + \varepsilon$), cálculo de inversas generalizadas $(X^TX)^-$, verificación de estimabilidad ($\lambda^T H = \lambda^T$) y resolución de ecuaciones normales.
2. **Formas Cuadráticas y ANOVA:** Obtención matricial de sumas de cuadrados ($SCReg$, $SCE$, $SCTotal$), cálculo del estimador insesgado de la varianza residual ($\hat{\sigma}^2$) y aplicación de propiedades de formas cuadráticas bajo normalidad.
3. **Medias Ponderadas e Intervalos de Confianza:** Estimación puntual y por intervalos ($1-\alpha$) para combinaciones lineales de tratamientos ponderadas por pesos poblacionales.
4. **Contrastes Lineales y Ortogonalidad:** Construcción de contrastes ortogonales independientes bajo la métrica definida por la matriz $(X^TX)^-$.
5. **Pruebas de Hipótesis Lineales Generales ($H_0: \Lambda\beta = d$):** Verificación de estimabilidad de $\Lambda$, cálculo de la suma de cuadrados de la hipótesis ($SCH_0$) y contraste de hipótesis conjuntas mediante el estadístico $F$.

---

## Estructura del Repositorio

```text
LinearModels_GM_Infe_applications/
├── Draft/
│   └── MLI.R                 # Script en R con los prototipos de cálculos matriciales
├── Informe/
│   ├── informe_final.Rmd     # Documento reproducible en R Markdown
│   ├── informe_final.pdf     # Informe final compilado con portada y tablas
│   ├── preamble.tex          # Configuración tipográfica y paquetes de LaTeX
│   └── logo.png              # Escudo institucional UNALM
├── Presentacion/
│   └── presentacion.html     # Material de exposición complementario
└── README.md                 # Documentación del proyecto
```
---

## Información Académica y Autores

- **Institución:** Universidad Nacional Agraria La Molina (UNALM)
- **Facultad:** Economía y Planificación (FEP) — Departamento Académico de Estadística e Informática
- **Curso:** Modelos Lineales I
- **Docente:** Mg. Samuel Huamaní Flores
- **Grupo 5:**
  - Jonnathan Jesús Pedraza Laboriano (20231505)
  - Luis Cortegana Quispe (20230388)
  - Miguel Encinas Alvarez (20230391)
  - Fiorella Fuente Bueno (20230395)
  - name (IDcode)