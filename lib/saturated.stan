data {
  int N;
  int SP;
  int Y[N, SP];
}

parameters {
  matrix[N, SP] W;
}

model {
  to_vector( W ) ~ normal(0.0,10.0);
  for (n in 1:N) {
    Y[n,] ~ multinomial( softmax(to_vector(W[n])) );
  }
}