Overview
========

This repository contains two adaptive numerical algorithms which use the Ogata quadrature method for applications to CSS TMD physics. The adaptive integration is optimized for integrals of the following form:


.. math::

  \begin{align}
  W_\nu(q_\perp;\vec{\alpha}) = \int_0^{\infty} \frac{db_\perp}{2\pi} b_\perp^{\nu+1} \widetilde{W}(b_\perp,\vec{\alpha}) J_\nu(b_\perp q_\perp).
  \end{align}

Where :math:`\vec{\alpha}` contains external kineamtic parameters and :math:`\widetilde{W}(b_\perp,\vec{\alpha})` is taken to be analytic close the the real axis. 

We denote the first Ogata quadrature method the 'untransformed adaptive Ogata' (agogu) and the second the 'transformed adaptive Ogata' (agogt). The untransformed Ogata quadrature is given by

.. math::

  \begin{align}
  \int_0^{\infty} {db_\perp}{2\pi} b_\perp^{\nu+1} \widetilde{W}(b_\perp,\vec{\alpha}) J_\nu(b_\perp q_\perp) \approx \frac{h}{2\pi}\frac{1}{q^{\nu+2}}\sum_{k = 1}^{N}\omega_{\nu k}x_{\nu k}^{\nu+1} \widetilde{W}(\frac{b_\perp}{q_\perp},\vec{\alpha}) J_\nu(x_{\nu k}).
  \end{align} 

The transformed quadrature sum is given by

.. math::

  \begin{align}
  \int_0^{\infty} {db_\perp}{2\pi} b_\perp^{\nu+1} \widetilde{W}(b_\perp,\vec{\alpha}) J_\nu(b_\perp q_\perp) \approx \frac{1}{2}\frac{1}{q^{\nu+2}}\sum_{k = 1}^{N} \omega_{\nu k} \psi'(x_{\nu k}) (\frac{\pi}{h}\psi(x_{\nu k}))^{\nu+1} \widetilde{W}(\frac{\pi}{h q_\perp}\psi(x_{\nu k}),\vec{\alpha}) J_\nu(\frac{\pi}{h}\psi(x_{\nu k}))
  \end{align}

Here 

.. math::

  \begin{align}
  w_{\nu k} = \frac{2}{\pi^2 \xi_{\nu |k|}J_{\nu+1}(\pi\xi_{\nu |k|})}\;,
  \qquad
  x_{\nu k} = h \xi_{\nu k}\; ,
  \qquad
  J_{\nu}(\pi \xi_{\nu k}) = 0\; .
  \end{align}

and :math:`h` and :math:`N` are the optimized parameters. See directed publication for specific details.
