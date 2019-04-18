Overview
========

This repository contains two adaptive numerical algorithms, written in python 2.7, which use the Ogata quadrature method for applications to CSS TMD physics. The adaptive integration is optimized for integrals of the following form


.. math::

  \begin{align}
  W_\nu(q_\perp;\vec{\alpha}) = \int_0^{\infty} \frac{db_\perp}{2\pi} b_\perp^{\nu+1} \widetilde{W}(b_\perp,\vec{\alpha}) J_\nu(b_\perp q_\perp).
  \end{align}

:math:`\vec{\alpha}` contains the external parameters and :math:`\widetilde{W}(b_\perp,\vec{\alpha})` is taken to be analytic close the the real axis and have a single peak at :math:`b_\perp>0`. The first algorithm uses equations 7 in the reference paper and is called the 'untransformed adaptive Ogata' method (adogu). The untransformed quadrature sum is given by

.. math::

  \begin{align}
  W_\nu(q_\perp;\vec{\alpha}) \approx \frac{h}{2\pi}\frac{1}{q^{\nu+2}}\sum_{k = 1}^{N}\omega_{\nu k}x_{\nu k}^{\nu+1} \widetilde{W}(\frac{x_{\nu k}}{q_\perp},\vec{\alpha}) J_\nu(x_{\nu k}).
  \end{align} 

The second algorithm uses equations 8 in the reference paper and is called the 'transformed adaptive Ogata' method (adogt). The transformed quadrature sum is given by

.. math::

  \begin{align}
  W_\nu(q_\perp;\vec{\alpha}) \approx \frac{1}{2}\frac{1}{q^{\nu+2}}\sum_{k = 1}^{N} \omega_{\nu k} \psi'(x_{\nu k}) (\frac{\pi}{h}\psi(x_{\nu k}))^{\nu+1} \widetilde{W}(\frac{\pi}{h q_\perp}\psi(x_{\nu k}),\vec{\alpha}) J_\nu(\frac{\pi}{h}\psi(x_{\nu k}))
  \end{align}

The parameters :math:`h` is optimized for adogu and adogt using equations 11 and 13 in the reference, respectively, and :math:`\psi(x) = x\tanh\left[\frac{\pi}{2}\sinh(x)\right]`. Here

.. math::

  \begin{align}
  w_{\nu k} = \frac{2}{\pi^2 \xi_{\nu |k|}J_{\nu+1}(\pi\xi_{\nu |k|})}\;,
  \qquad
  x_{\nu k} = h \xi_{\nu k}\; ,
  \qquad
  J_{\nu}(\pi \xi_{\nu k}) = 0\; .
  \end{align}

The integrator takes as an input the parameters :math:`\nu`, the order of the Bessel function, :math:`q_\perp`:, the transverse momentum, :math:`Q`, the estimated value of the reciprocal of the location of the peak of the function in :math:`b_\perp` space, and :math:`N`, the number of nodes. The function returns a single numerical value which gives an estimation of the desired integral. If :math:`Q` is not given, it is taken to be 10.
