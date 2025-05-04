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

%% Implementation of the scalar version using matrix arithmetic for the filtering and state updates
%% in order to remove the data dependencies.
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

%% Same as the matrix implementation, but instead of computing all outputs at once we pre-compute
%% the filter coefficients and unroll the result to compute two outputs at once.
function y = multival (x, Fn, Q)
  y = zeros(1, length(x));

  %% Compute feedforward/feedback coefficients
  fb = 1 / Q
  ff = tan(Fn * pi / 2)

  %% Compute overall gain
  gain = 1 / (1 + ff / Q + ff * ff)

  %% Compute the coefficient vectors.
  hpf = gain .* [1, -(ff + fb), -1, 0 ];
  bpf = gain .* [ff, 1, -ff, 0 ];
  lpf = gain .* [ff * ff, ff, 1 + fb, 0 ];

  %% Compute the state update vectors.
  s1 = ff .* hpf + bpf;
  s2 = ff .* bpf + lpf;

  %% Unroll.
  a = bpf;  % TODO: other filte rshapes
  b =  a(2) .* s1 +  a(3) .* s2 + [0, 0, 0,  a(1)];
  c = s1(2) .* s1 + s1(3) .* s2 + [0, 0, 0, s1(1)];
  d = s2(2) .* s1 + s2(3) .* s2 + [0, 0, 0, s2(1)];

  %% Do the filtering
  s = [0, 0]
  for m = 0:(length(x)/2 - 1)
    n = m * 2 + 1;
    z = [x(n), s(1), s(2), x(n + 1)];
    y(n)     = dot(a, z);
    y(n + 1) = dot(b, z);
    s(1)     = dot(c, z);
    s(2)     = dot(d, z);
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
y = multival(x, 0.5, 1.0);

% plot hp/bp/lp outputs
Y  = magnitude(y);
f  = linspace(0,1,length(Y));
plot(f, Y);
axis([0, 1, -24, 24])
