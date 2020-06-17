Overview
========

The algorithm in this repository optimizes Ogata's quadrature formula for TMD phenomenology purposes. The work in this repository is based on the Ogata quadrature formula which is given by

.. math::

  \begin{align}
  \int_{0}^{\infty} dx f(x) J_n(x) \approx \pi \sum_{j = 1}^{N} \omega_{n j}  f\left(\frac{\pi}{h}\psi(x_{n j})\right) J_n\left(\frac{\pi}{h}\psi(x_{n j})\right) \psi'(x_{n j})\,.
  \end{align}

:math:`w_{n j}` and :math:`x_{n j}` are the weights and nodes of the quadrature which are defined by

.. math::

  \begin{align}
  w_{n j} = \frac{2}{\pi^2 \xi_{n j}J_{n+1}(\pi\xi_{n j})}\;,
  \qquad
  x_{n j} = h \xi_{n j}\; ,
  \qquad
  J_n(\pi \xi_{n j}) = 0\;
  \end{align}

while the function :math:`\psi(t) = t \tanh\left(\frac{\pi}{2}\sinh(t)\right)`. The function :math:`J_n(x)` is the Bessel function of the first kind and :math:`n` is the order of the Bessel function. Here the function :math:`f(x)` must be even with respect to the y-axis. :math:`N` is the number of nodes used in the quadrature while :math:`h` is the node spacing parameter. For further details on this quadrature formula, please see reference [54] in our paper.

For TMD phenomenology, numerical Hankel transforms must be performed from from :math:`b_\perp`-space to transverse momentum space, :math:`q_\perp`-space. However, integrals of the form

.. math::

  \begin{align}
  W_\nu(q_\perp) = \int_0^{\infty} \frac{db_\perp}{2\pi} b_\perp^{\nu+1} \widetilde{W}(b_\perp) J_\nu(b_\perp q_\perp).
  \end{align}

tend to be huge bottlenecks of numerical computations.

The algorithm that we present exploits the fact that for TMD applications, the function :math:`b_\perp^{\nu+1} \widetilde{W}(b_\perp)` is a unimodal function. Our integrator takes takes as input :math:`\widetilde{W}(b_\perp)`, :math:`q_\perp`, :math:`N`, as well as :math:`Q`, a guess for the inverse of the value of :math:`b_\perp` at which the function :math:`b_\perp^{\nu+1} \widetilde{W}(b_\perp)` is maximized. From this function and these parameters, the algorithm uses Eq. (21) and Eq.(23) in our paper to find that value of :math:`h` which optimizes the performance of the quadrature.
