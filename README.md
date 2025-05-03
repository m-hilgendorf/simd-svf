# simd-svf

Experimenting with simd optimizations in a state-variable filter.

## Derivation
### Tl;Dr

A ZDF/TPT state variable filter can be computed as a pair of matrix multiplications, which removes all data dependencies when evaluating the filter:

#### Inputs
$$
  F_n: \text{ Normalized Frequency} \\
  Q: \text{ Q factor } \\
  s_1, s_2: \text{state variables, initialized to 0}
$$
#### Formula
$$
  \text{ff} = \tan(F_n \frac{pi}{2})
$$
$$
  \text{fb} = \frac{1}{Q}
$$
$$
  \text{G}  = \frac{1}{1 + (\text{ff } \text{fb}) + \text{ff}^2}
$$
$$
  A = G \begin{pmatrix}
  1 & -(\text{ff } + \text{fb }) & -1 \\
  \text{ff } & 1 & -\text{ff } \\
  \text{ff }^2 & \text{ff } & (1 + \text{fb})
  \end{pmatrix}
$$
$$
B = \begin{pmatrix}
  \text{ff } A_{11...13} + A_{21...23} \\
  \text{ff } A_{21...23} + A_{31...33}
\end{pmatrix}
$$
$$
\begin{pmatrix}
y_{hp} \\
y_{bp} \\
y_{lp} \\
\end{pmatrix}
= A
\begin{pmatrix}
  x \\
  s_1 \\
  s_2
\end{pmatrix}
$$
$$
\begin{pmatrix}
s_1 \\
s_2
\end{pmatrix}
= B
\begin{pmatrix}
  x \\
  s_1 \\
  s_2
\end{pmatrix}
$$

### Control inputs
In [the Art of VA Filter Design](https://www.native-instruments.com/fileadmin/ni_media/downloads/pdf/VAFilterDesign_2.1.0.pdf), Vadim Zavalishin gives a discretization of a 2nd order SVF. He uses some different variable namings, so I will use different terms. Our control inputs are normalized frequency and Q factor.

$$
F_n : \text{Normalized Frequency} = \omega_c T/ 2\pi \\
Q: \text{Q factor} = \frac{1}{2R}
$$

I prefer normalized frequency to radial frequency and sampling interval as it makes the internals sample-rate agnostic. I prefer Q factor because it is more common on EQs that I have used than an arbitrary "resonance" parameter that you may find on a synthesizer.

## Internal coefficients

There are two internal coefficients worth computing up front, which I will call $\text{ff}$ and $\text{fb}$ for "feed forward" and "feed back". The feed forward coefficient is the gain applied at each integrator (what Zavalishin refers to as $g$) and the feedback coefficient is the gain applied to the feedback path from the bandpass output ($2R = 1/Q$).

$$
  \text{ff} = \tan(F_n \frac{\pi}{2}) \\
  \text{fb} = \frac{1}{Q}
$$

### Base derivation

Zavalishin derives the system of equations for the highpass, bandpass, and lowpass outputs of the SVF:

$$
y_{hp} = x - \text{fb } y_{bp} - y_{lp} \\
y_{bp} = \text{ff } y_{hp} + s_1 \\
y_{lp} = \text{ff } y_{bp} + s_2
$$

The problem with this system of equations is that it contains zero-delay feedback loops. To realize it in code we need to solve this system of equations for one of the $y$ variables.

In the book, there are derivations in terms of $y_{hp}$ and $y_{bp}$. We go ahead and solve for all three.

$$
y_{hp} = \frac
{x - (\text{ff} + \text{fb}) s_1 - s_2 }
{1 + \text{ff } \text{fb } + \text{ff}^2}
$$

$$
y_{bp} = \frac
{\text{ff } x + s_1 - \text{ff }s_2 }
{1 + \text{ff } \text{fb } + \text{ff}^2}
$$

$$
y_{lp} = \frac
{\text{ff }^2 x + \text{ff } s_1 + (1 + \text{fb}) s_2 }
{1 + \text{ff } \text{fb } + \text{ff}^2}
$$

We can pull out the shared constant and call it $G = \frac{1} {1 + \text{ff } \text{fb } + \text{ff}^2}$

We can also state this as a matrix multiplication

$$
\begin{pmatrix}
y_{hp} \\
y_{bp} \\
y_{lp} \\
\end{pmatrix}
= G
\begin{pmatrix}
1 & -(\text{ff } + \text{fb }) & -1 \\
\text{ff } & 1 & -\text{ff } \\
\text{ff }^2 & \text{ff } & (1 + \text{fb})
\end{pmatrix}
\begin{pmatrix}
x \\
s_1 \\
s_2
\end{pmatrix}
$$

The matrix formulation is useful because it removes all data dependencies when computing the outputs, which on paper should improve instruction-level parallelism.

We can also formulate the state update as a matrix multiplication using the rows of the coefficient matrix above, call it $A$.

$$
\begin{pmatrix}
  s_1 \\
  s_2 \\
\end{pmatrix}
=
\begin{pmatrix}
  \text{ff } A_{11...13} + A_{21...23} \\
  \text{ff } A_{21...23} + A_{31...33}
\end{pmatrix}
\begin{pmatrix}
  x \\
  s_1 \\
  s_2 \\
\end{pmatrix}
$$
