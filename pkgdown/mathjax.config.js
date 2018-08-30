window.MathJax = {
  tex2jax: {
    inlineMath:  [ ['\\(','\\)'], ['$', '$'] ],
    displayMath: [ ['\\[','\\]'], ['$$', '$$'] ],
    processEscapes: true
  },
  TeX: {
    Macros: {
      a: "\\mathbf{a}",
      A: "\\mathbf{A}",
      B: "\\mathbf{B}",
      C: "\\mathbf{B}",
      D: "\\mathbf{D}",
      e: "\\mathbf{e}",
      E: "\\mathbf{E}",
      F: "\\mathbf{F}",
      g: "\\mathbf{g}",
      H: "\\mathbf{H}",
      I: "\\mathbf{I}",
      L: "\\mathbf{L}",
      M: "\\mathbf{M}",
      P: "\\mathbf{P}",
      Q: "\\mathbf{Q}",
      R: "\\mathbf{R}",
      u: "\\mathbf{u}",
      U: "\\mathbf{U}",
      v: "\\mathbf{v}",
      V: "\\mathbf{V}",
      w: "\\mathbf{w}",
      W: "\\mathbf{W}",
      x: "\\mathbf{x}",
      X: "\\mathbf{X}",
      y: "\\mathbf{y}",
      Y: "\\mathbf{Y}",
      z: "\\mathbf{z}",
      Z: "\\mathbf{Z}",
      one: "\\mathbf{1}",
      zero: "\\mathbf{0}",

      xx: "\\tilde{\\mathbf{x}}",
      XX: "\\tilde{\\mathbf{X}}",
      yy: "\\tilde{\\mathbf{y}}",
      YY: "\\tilde{\\mathbf{Y}}",
      zz: "\\tilde{\\mathbf{z}}",
      ZZ: "\\tilde{\\mathbf{Z}}",
      
      cR: "\\mathcal{R}",
      cN: "\\mathcal{N}",

      bb: "\\boldsymbol{\\beta}",
      bh: "\\hat{\\beta}",
      bbh: "\\boldsymbol{\\hat{\\beta}}",
      gam: "\\gamma",
      eps: "\\epsilon",
      veps: "\\varepsilon",
      be: "\\boldsymbol{\\eta}",
      lam: "\\lambda",
      bL: "\\boldsymbol{\\Lambda}",
      bm: "\\boldsymbol{\\mu}",
      bS: "\\boldsymbol{\\Sigma}",
      th: "\\hat{\\theta}",
      bt: "\\boldsymbol{\\theta}",
      bth: "\\boldsymbol{\\hat{\\theta}}",
   
      Pr: "\\mathbb{P}",
      Ex: "\\mathbb{E}",
      Var: "\\mathbb{V}",
      cov: "\\textrm{Cov}",
      Norm: "\\textrm{N}", 
      Pois: "\\textrm{Pois}", 
      Exp: "\\textrm{Exp}",
      trace: "\\textrm{tr}",
      ind: "\\perp\\!\\!\\!\\perp",
   
      inD: "\\overset{\\text{d}}{\\longrightarrow}", 
      inP: "\\overset{\\text{p}}{\\longrightarrow}", 
      adist: "\\overset{\\text{.}}{\\sim}",
      SD: "\\textrm{SD}",
      KL: "\\textrm{KL}",

      inprod: ['\\langle #1, #2 \\rangle', 2],

      norm: ["\\left\\lVert #1 \\right\\rVert", 1],
      abs: ["\\left\\lvert #1 \\right\\rvert", 1],
      as: ['\\begin{align*} #1 \\end{align*}', 1],
      als: ['\\begin{align}\\label{#1}\\begin{split}#2\\end{split}\\end{align}', 2],
      al: ['\\begin{align}\\label{#1} #2 \\end{align}', 2]
    },
    equationNumbers: { autoNumber: "AMS" },
    extensions: ["autoload-all.js"]
  }
};
