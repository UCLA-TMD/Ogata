Overview
========

This repository contains the optimized Ogata quadrature numerical algorithm, available in both python 2.7 and Fortran 77. This algorithm uses the Ogata quadrature method for applications to CSS TMD physics. The integration optimizes integrals of the following form


.. math::

  \begin{align}
  W_\nu(q_\perp) = \int_0^{\infty} \frac{db_\perp}{2\pi} b_\perp^{\nu+1} \widetilde{W}(b_\perp) J_\nu(b_\perp q_\perp).
  \end{align}

Here  :math:`\widetilde{W}(b_\perp)` is a function which is taken be  analytic close the the real axis and have a single peak at :math:`b_\perp>0`. The algorithm uses equations 15 in the reference paper and is called the 'optimized Ogata' method. The quadrature sum is given by

.. math::

  \begin{align}
  W_\nu(q_\perp) \approx \frac{1}{2q_\perp}\sum_{k = 1}^{N} \omega_{\nu k} \psi'(x_{\nu k}) (\frac{\pi}{h q_\perp}\psi(x_{\nu k}))^{\nu+1} \widetilde{W}(\frac{\pi}{h q_\perp}\psi(x_{\nu k})) J_\nu(\frac{\pi}{h}\psi(x_{\nu k}))
  \end{align}

The parameters :math:`h` is optimized using equations (19) and (21) in the reference. Furthermore :math:`\psi(x) = x\tanh\left[\frac{\pi}{2}\sinh(x)\right]`. Here

.. math::

  \begin{align}
  w_{\nu k} = \frac{2}{\pi^2 \xi_{\nu |k|}J_{\nu+1}(\pi\xi_{\nu |k|})}\;,
  \qquad
  x_{\nu k} = h \xi_{\nu k}\; ,
  \qquad
  J_{\nu}(\pi \xi_{\nu k}) = 0\; .
  \end{align}

The integrator takes as an input the parameters :math:`\nu`, the order of the Bessel function, :math:`q_\perp`:, the transverse momentum, :math:`Q`, the estimated value of the reciprocal of the location of the peak of the function in :math:`b_\perp` space, and :math:`N`, the number of nodes. The function returns a single numerical value which gives an estimation of the desired integral. If :math:`Q` is not given, it is taken to be 10.
