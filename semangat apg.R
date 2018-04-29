###############################################################################
####################### Rata2, Varians, Kovarians, Korelasi ###################
###############################################################################
mu = function(data){
  data = as.array(data)
  n = length(data)
  y = 0
  for (item in data){
    y = y + item
  }
  return(y/n)
}

mean_vector = function(data){
  mat = matrix(nrow = dim(data)[2],ncol = 1)
  for(i in 1:dim(data)[2]){
    mat[i,1] = mu(data[,i])
  }
  return(mat)
}

varians = function(data){
  n = length(data)
  y = 0
  m = mu(data)
  for (item in data){
    y = y + (item-m)^2
  }
  return(y/n)
}

covarians = function(data1,data2){
  if (length(data1) != length(data2)){
    print("Panjang data tidak sama, tidak bisa menghitung kovarians")
  } else {
    n = length(data1)
    y = 0
    m1 = mu(data1)
    m2 = mu(data2)
    for (i in 1:n){
      y = y + (data1[i]-m1)*(data2[i]-m2)
    }
    return(y/(n-1))
  }
}

cov_matrix = function(data){
  mat = matrix(nrow = dim(data)[2],ncol = dim(data)[2])
  for(i in 1:dim(data)[2]){
    for(j in 1:dim(data)[2]){
      mat[i,j] = covarians(data[,i], data[,j])
    }
  }
  return(mat)
}

correlation = function(data1,data2){
  rho = covarians(data1,data2)/sqrt(varians(data1)*varians(data2))
  return(rho)
}

corr_matrix = function(data){
  mat = matrix(nrow = dim(data)[2],ncol = dim(data)[2])
  for(i in 1:dim(data)[2]){
    for(j in 1:dim(data)[2]){
      mat[i,j] = correlation(data[,i], data[,j])
    }
  }
  return(mat)
}

###############################################################################
############################## Normal Distribution ############################
###############################################################################

qq_plot_uni = function(data){
  data = sort(data)
  n = length(data)
  Pz = c()
  for( i in 1:n){
    z = (i-0.5)/n
    Pz = c(Pz,pnorm(z))
  }
  par(pty="s")
  plot(Pz,data)
}

d_2 = function(data){
  n = dim(data)[1]
  m = mean_vector(data)
  S = cov_matrix(data)
  Dj = c()
  for (i in 1:n){
    xj = t(t(data[i,]))
    dj = t(xj-m)%*%solve(S)%*%(xj-m)
    Dj = c(Dj, dj)
  }
  return(Dj)
}

qq_plot_multi = function(data){
  v = dim(data)[2]
  Dj = sort(d_2(data))
  n = length(Dj)
  Pchi = c()
  for (i in 1:n){
    chi = (i-0.5)/n
    Pchi = c(Pchi, pchisq(chi,v))
  }
  par(pty="s")
  plot(Pchi, Dj)
}

################################################################################
################### Pengujian Hipotesis Beda Rata - Rata #######################
################################################################################

uji_vektor_rata2_cov = function(data, mcov, mmu, alpha){
  m = mean_vector(data)
  n = dim(data)[1]
  v = dim(data)[2]
  z2 = t(m-mmu)%*%solve(mcov/n)%*%(m-mmu)
  pval = 2*pt(-abs(z2), df = v)
  print("P-Value : ", pval)
  if (pval<=alpha){
    print("Keputusan Tolak H0")
  } else {
    print("Keputusan Gagal Tolak H0")
  }
}

uji_vektor_rata2_hotteling = function(data,mmu,alpha){
  m = mean_vector(data)
  n = dim(data)[1]
  
  T2 = t(m-mmu)%*%solve(cov_matrix(data)/n)%*%(m-mmu)
  alpha = 1-alpha
  p = length(m)
  v = n-1
  Fval = qf(alpha, p, v-p+1)
  if (T2 > (v*p/(v-p+1)*Fval)){
    print("Keputusan Tolak H0")
    print("Selang Kepercayaan : ")
    for(i in 1:p){
      y = sqrt(p*v/(n-p)*Fval*varians(data[,i])/n)
      b_bawah = m[i] - y
      b_atas = m[i] + y
      cat("CI-", i, " : ", b_bawah, "<= ?? <=", b_atas )
    }
  } else {
    print("Keputusan Gagal Tolak H0")
  }
}

perbandingan_vektor_rata2 = function(data1, data2, alpha){
  m1 = mean_vector(data1)
  m2 = mean_vector(data2)
  n1 = dim(data1)[1]
  n2 = dim(data2)[1]
  W1 = (n1-1)*cov_matrix(data1)
  W2 = (n2-1)*cov_matrix(data2)
  Sgab= (W1+W2)/(n1+n2-2)
  T2 = n1*n2/(n1+n2)%*%t(m1-m2)%*%solve(Sgab)%*%(m1-m2)
  cat("T2 : ", T2, "\n")
  p = length(m1)
  Fval = qf(alpha, p, n1+n2-p-1)
  c2 = (n1+n2-2)*p/(n1+n2-p-1)*Fval
  if (T2 > c2){
    print("Keputusan Tolak H0")
  } else {
    print("Keputusan Gagal Tolak H0")
  }
  print("CI Hotelling : ")
  for (i in 1:p){
    m = m1[i] - m2[i]
    b_bawah = m - c2*sqrt((n1+n2)/(n1*n2)*Sgab[i,i])
    b_atas = m + c2*sqrt((n1+n2)/(n1*n2)*Sgab[i,i])
    cat("CI-", i, " : ", b_bawah, "<= ??1-??2 <=", b_atas ,"\n" )
  }
  print("CI Bonferoni : ")
  tval = qt(alpha/(2*p), n1+n2-2)
  for(i in 1:p){
    m = m1[i] - m2[i]
    b_bawah = m - tval*sqrt((n1+n2)/(n1*n2)*Sgab[i,i])
    b_atas = m + tval*sqrt((n1+n2)/(n1*n2)*Sgab[i,i])
    cat("CI-", i, " : ", b_bawah, "<= ??1-??2 <=", b_atas ,"\n")
  }
}

uji_data_berpasangan = function(data1, data2){
  if (dim(data1)[1]!=dim(data2)[1]){
    print("Jumlah data tidak sama")
  } else {
    p = dim(data1)[2]
    n = dim(data1)[1]
    Di = c()
    for (i in 1:p){
      Di= cbind(Di, data1[,i]-data2[,i])
    }
    Di = as.matrix(Di)
    d_bar = mean_vector(Di)
    Sd = cov_matrix(Di)
    T2 = n*t(d_bar)%*%solve(Sd)%*%d_bar
    Fval = qf(alpha, 2, n-2)
    if (T2 > 2*(n-1)*Fval/(n-2)){
      print("Keputusan Tolak H0")
    } else {
      print("Keputusan Gagal Tolak H0")
    }
  }
}

################################################################################
############################# Profile Analysis #################################
################################################################################

mat_contrass = function(p){
  for (i in 1:(p-1)){
    row = c(rep(0,i-1), -1, 1)
    row = c(row, rep(0,(p-length(row))))
    print(row)
  }
}

PA_satu_sampel = function(data, alpha){
  
}
paralel_PA_dua_sampel = function(data1, data2, alpha){
  m1 = mean_vector(data1)
  m2 = mean_vector(data2)
  n1 = dim(data1)[1]
  n2 = dim(data2)[1]
  W1 = (n1-1)*cov_matrix(data1)
  W2 = (n2-1)*cov_matrix(data2)
  Sgab= (W1+W2)/(n1+n2-2)
  p = dim(data1)[2]
  C = mat_contrass(p)
  T2 = n1*n2/(n1+n2)*t(m1-m2)%*%t(C)%*%solve(C%*%Sgab%*%t(C))%*%C%*%(m1-m2)
  
}
mat_contrass(4)

################################################################################
############################# MANOVA ###########################################
################################################################################

manova_2_sampel = function(data1, data2){
  install.packages("knitr")
}