Overview
========

The algorithm in this repository optimizes Ogata's quadrature formula for TMD phenomenology. This work is based on the Ogata quadrature formula, which is given by the expression

.. math::

  \begin{align}
  \int_{0}^{\infty} dx f(x) J_n(x) \approx \pi \sum_{j = 1}^{N} \omega_{n j}  f\left(\frac{\pi}{h}\psi(x_{n j})\right) J_n\left(\frac{\pi}{h}\psi(x_{n j})\right) \psi'(x_{n j})\,.
  \end{align}

In this expression, :math:`w_{n j}` and :math:`x_{n j}` are the weights and nodes of the quadrature which are defined by

.. math::

  \begin{align}
  w_{n j} = \frac{2}{\pi^2 \xi_{n j}J_{n+1}(\pi\xi_{n j})}\;,
  \qquad
  x_{n j} = h \xi_{n j}\; ,
  \qquad
  J_n(\pi \xi_{n j}) = 0\;
  \end{align}

while 

.. math::
  \psi(t) = t \tanh\left(\frac{\pi}{2}\sinh(t)\right)

and the function :math:`J_n(x)` is the Bessel function of the first kind of order :math:`n`. The quadrature has two parameters which control the numerical errors, :math:`N`, the number of nodes used in the quadrature, and :math:`h`, is the node spacing parameter. For further details on this quadrature formula, please see reference [54] in our paper.

For TMD phenomenology, numerical Hankel transformations must be performed from :math:`b_\perp`-space to transverse momentum space, :math:`q_\perp`-space. However, integrals of the form

.. math::

  \begin{align}
  W(q_\perp) &= \int_0^{\infty} \frac{db_\perp}{2 \pi} b_\perp^{n+1} \widetilde{W}(b_\perp) J_n(b_\perp q_\perp) \\
                 &= \frac{1}{q_\perp^{n+2}}\int_0^{\infty} \frac{dx}{2 \pi} x^{n+1} \widetilde{W}\left(\frac{x}{q_\perp}\right) J_n(x)
  \end{align}

tend to be huge bottlenecks of numerical computations. Note that to go from the first to the second line, we have made the substitution :math:`x = b_\perp q_\perp` in order to match the form of the Ogata quadrature formula.

In order to optimize the quadrature formula, we exploit the fact that the function :math:`b_\perp^{n+1} \widetilde{W}(b_\perp)` is unimodal on the positive :math:`b_\perp`-axis. Our integrator takes as input :math:`\widetilde{W}(b_\perp)`, :math:`q_\perp`, :math:`N`, :math:`n`, as well as :math:`Q`, an initial guess for the inverse of the value of :math:`b_\perp` at which the function :math:`b_\perp^{n+1} \widetilde{W}(b_\perp)` is maximized. From this input, the algorithm uses Eq. (21) and Eq.(23) in our paper to find that value of :math:`h` which optimizes the performance of the quadrature. The integrator outputs the numerical value

.. math::

  \begin{align}
  \textrm{FBT}(n,\widetilde{W}\left(b_\perp\right),q_\perp,Q,N) \approx \frac{1}{2 q_\perp^{2n+2}} \sum_{j = 1}^{N} \omega_{n j} & \left(\frac{\pi}{h}\psi(x_{n j})\right)^{2n+1} \widetilde{W}\left(\frac{1}{q_\perp}\frac{\pi}{h}\psi(x_{n j})\right) \\ & \times J_n\left(\frac{\pi}{h}\psi(x_{n j})\right) \psi'(x_{n j})\,.
  \end{align}

It is important to note when using the integrator, the factor of :math:`\frac{1}{2\pi}` is included in the numerical output and that the function :math:`\widetilde{W}` must be :math:`b_\perp` dependent.
