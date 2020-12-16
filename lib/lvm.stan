data {
  int N;
  int E;
  int SP;
  int L;
  int Y[N, SP];
  matrix[N, E] X;
}

parameters {
  matrix[E, SP] W;
  matrix[N, L] LV;
  matrix[L, SP] LF;
}

model {
  matrix[N, SP] pred;
  pred = X * W + LV * LF;
  to_vector(W) ~ normal(0.0, 3.0);
  to_vector(LV) ~ normal(0.0, 3.0);
  to_vector(LF) ~ normal(0.0, 3.0);
  for (n in 1:N) {
    Y[n,] ~ multinomial( softmax(to_vector(pred[n])) );
  }
}