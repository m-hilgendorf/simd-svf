clear all

%% Reference scalar implementation
%%
%% Arguments:
%%  - x: input signal
%%  - Fn: normalized frequency.
%%  - Q: Q-factor (resonance = 1/2Q). Should clamped in the range [0.5; 10)
%%
function y = scalar (x, Fn, Q)
  %% Initialize y
  y = zeros(3, length(x));

  %% Compute feedforward/feedback coefficients
  fb = 1 / Q
  ff = tan(Fn * pi / 2)

  %% Compute overall gain
  gain = 1 / (1 + ff / Q + ff * ff)

  %% Initialize the state.
  s = [0, 0];

  %% Do the filtering.
  for n = 1:length(x)
    %% As per the Zavalishin text, solve in terms of yhp and then compute bp/lp outputs.
    yhp = gain * (x(n) - fb * s(1) - ff * s(1) - s(2))
    ybp = ff * yhp + s(1);
    ylp = ff * ybp + s(2);

    % Update the state variables.
    s(1) = ff * yhp + ybp;
    s(2) = ff * ybp + ylp;

    y(:,n) = [yhp, ybp, ylp];
    %% Update the st
  endfor
endfunction

%% Matrix implementation.
%%
%% Arguments:
%%  - x: input signal
%%  - Fn: normalized frequency.
%%  - Q: Q-factor (resonance = 1/2Q). Should clamped in the range [0.5; 10)
%%
%% Reference:
%%    This implementation is based on the discretization of the SVF found in "The Art of VA Filter
%%    Design" by Vadim Zavalishin. His TPT approach works by discretizing the integrators in the
%%    analog topology and then algebraically resolving the zero-delay-feedback loops by solving the
%%    system of equations by hand.
%%
%%    The matrix formulation is derived by solving for the highpass, bandpass, and lowpass outputs
%%    independently to create the coefficient matrix.
%%
%% There are a handful of useful optimizations that this prototype can be used to design:
%%  - The tan(Fn * pi / 2) call can be replaced with an approximation. Bad approximations show up
%%    as errors in the cutoff frequency, usually worse the closer to Nyquist.
%%  - There are no data dependencies on any variable (input, output, or state).
%%  - If the desired output is highpass/bandpass/lowpass, the matrix multiply Az can be replaced
%%    with a dot product of the corresponding shape dot(A(1), z) for highpass, A(2) for bandpass,
%%    A(3) for lowpass.
%%  - The matrix multiply can be "unrolled" to compute 2 outputs at the same time.
%%
function y = matrix (x, Fn, Q)
  %% Initialize y
  y = zeros(3, length(x));

  %% Compute feedforward/feedback coefficients
  fb = 1 / Q
  ff = tan(Fn * pi / 2)

  %% Compute overall gain
  gain = 1 / (1 + ff / Q + ff * ff)

  %% Initialize the state variables.
  s = [0, 0];

  %% Do the filtering.
  for n = 1:length(x)
    %% The "A" matrix maps input and state variables to highpass/lowpass/bandpass outputs. These
    %% matrices can be computed ahead of time if the filter is not time-variant.
    A = gain .* [
      1, -(ff + fb), -1 ;
      ff, 1, -ff ;
      ff * ff, ff, 1 + fb ;
    ];

    %% The "B" matrix maps state to next state, and it is derived from the "A" matrix.
    B = [
      ff .* A(1, :) + A(2, :),
      ff .* A(2, :) + A(3, :),
    ];

    %% Create the input vector, by concatenating (x) with the state vector.
    z = [x(n), s(1), s(2)];

    %% Compute output and next state. Note that there is no data dependency between the two.
    y(:, n) = A * z';
    s       = B * z';
  endfor
endfunction

function H = magnitude(y)
  H = fft(y)(1:length(y)/2);
  H = 20*log10(abs(H));
endfunction

% unit implse
x = zeros(1, 512);
x(1) = 1;

% filter with some parameters
y = matrix(x, 0.5, 5.0);

% plot hp/bp/lp outputs
hp = magnitude(y(1, :));
bp = magnitude(y(2, :));
lp = magnitude(y(3, :));
f  = linspace(0,1,length(hp));
plot(f, hp, f, bp, f, lp);
axis([0, 1, -24, 24])
