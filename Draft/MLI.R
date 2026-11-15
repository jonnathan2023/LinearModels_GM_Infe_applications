
consumo <- c(74, 68, 77, 76, 80, 85, 93)
educacion <- c("1", "1", "1", "2", "2", "3", "3")

datos <- data.frame(educacion, consumo)
datos$educacion <- as.factor(datos$educacion)
datos

# Matriz de diseño

modelo <- lm(consumo ~ educacion - 1, data = datos)
X <- model.matrix(modelo)
X <- cbind(Intercepto = rep(1, length(consumo)), X)
colnames(X) <- NULL
X

y <- as.matrix(consumo); y

# Cálculo de la inversa generalizada

# Función para calcular inversa generalizada
invg <- function(X) {
  p <- ncol(A) # Número de columnas (dimensión de A)
  r <- qr(A)$rank # Rango de X'X
  n <- p - r # Número de filas/columnas dependientes
  m <- n + 1
  A_reducida <- A[-(1:n), -(1:n)] # Extraer submatriz no singular
  D_inv <- solve(A_reducida) # Inversa de la parte no singular
  G <- matrix(0, nrow = p, ncol = p) # Inicializar matriz G con ceros
  G[m:p, m:p] <- D_inv # Insertar la inversa en el bloque inferior derecho
  return(G)
}

A <- t(X) %*% X; A

library(MASS)
G <- invg(A); fractions(G)

# b. Encuentre una solución particular de la ecuación normal.

b0 <- G %*% t(X) %*% y; b0

# c. ANOVA

n <- length(y); n

SCReg <- t(y) %*% X %*% G %*% t(X) %*% y; SCReg

SCTotal <- t(y) %*% y; SCTotal

SCMedia <- n*mean(y)^2; SCMedia

SCResid <- SCTotal - SCReg; SCResid

SCReg_m <- SCReg - SCMedia; SCReg_m

s2 <- SCResid/(n - qr(X)$rank); s2

CMReg_m <- SCReg_m/(qr(X)$rank - 1); CMReg_m

Fcal <- CMReg_m/s2; Fcal

qf(0.95, qr(X)$rank - 1, n - qr(X)$rank)

# d. Es estimable \mu + \alpha? Y si es así, calcule su B.L.U.E (Media y
#  Varianza) e intervalo de confianza

l <- matrix(c(1, 1, 0, 0), nrow = 4); l

H <- G %*% t(X) %*% X; H

all.equal(t(l), t(l) %*% H)

# Media
t(l) %*% b0

# Varianza
s2 %*% (t(l) %*% G %*% l)

# Intervalo de confianza
qt(0.975, n - qr(X)$rank)

LI <- (t(l) %*% b0) - qt(0.975, n - 
                           qr(X)$rank) * sqrt(s2 %*% (t(l) %*% G %*% l))
LI

LS <- (t(l) %*% b0) + qt(0.975, n - 
                           qr(X)$rank) * sqrt(s2 %*% (t(l) %*% G %*% l))
LS

# e.

Lambda <- matrix(c(0, -1, 1, 0,
                   0, -1, -1, 2), nrow = 2, byrow = TRUE); Lambda

d <- matrix(c(9, 30), nrow = 2); d

# Ver si es estimable

all.equal(Lambda %*% H, Lambda)

qr(Lambda)$rank

SCH0 <- t(Lambda %*% b0 - d) %*% solve(Lambda %*% G %*% t(Lambda)) %*%
  (Lambda %*% b0 - d); SCH0

Fcal <- SCH0/(2 * s2); Fcal
Fcrit <- qf(0.95, 2, n - qr(X)$rank); Fcrit

# f.

Lambdaf <- matrix(c(7, 3, 2, 2), nrow = 1, byrow = TRUE); Lambdaf

# Ver si es estimable

all.equal(Lambdaf %*% H, Lambdaf)

qr(Lambdaf)$rank

SCH0_f <- t(Lambdaf %*% b0) %*% solve(Lambdaf %*% G %*% t(Lambdaf)) %*%
  (Lambdaf %*% b0); SCH0_f

Fcal_f <- SCH0_f/(1 * s2); Fcal_f
Fcrit_f <- qf(0.95, 1, n - qr(X)$rank); Fcrit_f

# g.

Lambdag <- matrix(c(0, 1, -1, 0,
                    0, 1, 0, -1), nrow = 2, byrow = TRUE); Lambdag

# Ver si es estimable

all.equal(Lambdag %*% H, Lambdag)

qr(Lambdag)$rank

SCH0_g <- t(Lambdag %*% b0) %*% solve(Lambdag %*% G %*% t(Lambdag)) %*%
  (Lambdag %*% b0); SCH0_g

Fcal_g <- SCH0_g/(2 * s2); Fcal_g
Fcrit_g <- qf(0.95, 2, n - qr(X)$rank); Fcrit_g

# h.

Lambdah1 <- matrix(c(7, 3, 2, 2), nrow = 4); Lambdah1
Lambdah2 <- matrix(c(0, 1, -1, 0), nrow = 4); Lambdah2
Lambdah3 <- matrix(c(0, 1, 0, -1), nrow = 4); Lambdah3

Lambdah <- matrix(c(Lambdah1, Lambdah2, Lambdah3), nrow = 3, byrow = TRUE)
Lambdah

## 1. Estimabilidad

all.equal(t(Lambdah1), t(Lambdah1) %*% H)
all.equal(t(Lambdah2), t(Lambdah2) %*% H)
all.equal(t(Lambdah3), t(Lambdah3) %*% H)

# Independencia Lineal

qr(Lambdah)$rank

# 2. Ortoganalidad

t(Lambdah1) %*% G %*% Lambdah2
t(Lambdah1) %*% G %*% Lambdah3
t(Lambdah2) %*% G %*% Lambdah1
t(Lambdah2) %*% G %*% Lambdah3
t(Lambdah3) %*% G %*% Lambdah1
t(Lambdah3) %*% G %*% Lambdah2

# 2.

## 0.12 <- Nivel 1 (\mu + \alpha_1)
## 0.68 <- Nivel 2 (\mu + \alpha_2)
## 0.20 <- Nivel 3 (\mu + \alpha_3)

## La media poblacional es el promedio ponderaso de las medias de
## cada grupo

## \phi = (0.12 + 0.68 + 0.20)\mu + 0.12\alpha_1 + 0.68\alpha_2 + 0.20\alpha_3
## \phi = \mu + 0.12\alpha_1 + 0.68\alpha_2 + 0.20\alpha_3

## \lambda^T = (1, 0.12, 0.68, 0.20)

## \lambda^T es estimable porque \lambda_1 + \lambda_2 + \lambda_3 = \lambda_0

# 2.a. Media poblacional estimadad (\lambda^T \hat{b})

lambda2_0 <- matrix(c(1, 0.12, 0.68, 0.20), nrow = 4); lambda2_0 

mediaestimada2 <- t(lambda2_0) %*% b0; mediaestimada2

# 2.b. Varianza estimada del estimador

## \hat{Var(\lambda^T \hat{\beta})} = \hat{\sigma}^2 \lambda^T G \lambda

varestimada2 <- s2 %*% t(lambda2_0) %*% G %*% lambda2_0; varestimada2

# 2.c. Intervalo de confianza simétrico al 95% para el promedio poblacional

n <- length(y); n
mediaestimada2 - qt(1 - 0.05/2, n - qr(X)$rank)*sqrt(varestimada2)
mediaestimada2 + qt(1 - 0.05/2, n - qr(X)$rank)*sqrt(varestimada2)

# 2.d.

Lambda_2d <- matrix(c(1, 1, 0, 0,
                      0, 1, 0, -1), nrow = 2, byrow = TRUE); Lambda_2d

d_2d <- matrix(c(70, -15), nrow = 2); d_2d

# Ver si es estimable

all.equal(Lambda_2d %*% H, Lambda_2d)

qr(Lambda_2d)$rank

SCH0_2d <- t(Lambda_2d %*% b0 - d_2d) %*% solve(Lambda_2d %*% G %*% t(Lambda_2d)) %*%
  (Lambda_2d %*% b0 - d_2d); SCH0_2d

Fcal_2d <- SCH0_2d/(2 * s2); Fcal_2d
Fcrit_2d <- qf(0.95, 2, n - qr(X)$rank); Fcrit_2d

# 2.e.



############################################
############################################
# 3.

puntaje <- c(37, 25, 42, 28,
             23, 29,
             38, 30, 25,
             23, 29)
trabajadores <- c(rep("Obrero", 4), rep("Artesanos", 2),
               rep("Profesionales", 3), rep("Independientes", 2))

datos <- data.frame(puntaje, trabajadores)
datos$trabajadores <- factor(datos$trabajadores,
                             levels = c("Obrero", "Artesanos",
                                        "Profesionales", "Independientes"))
datos

# Matriz de diseño

modelo <- lm(puntaje ~ trabajadores - 1, data = datos)
X3 <- model.matrix(modelo)
X3 <- cbind(Intercepto = rep(1, length(puntaje)), X3)
colnames(X3) <- NULL
X3

y3 <- as.matrix(puntaje); y3

# Cálculo de la inversa generalizada

# Función para calcular inversa generalizada
invg <- function(X) {
  p <- ncol(A) # Número de columnas (dimensión de A)
  r <- qr(A)$rank # Rango de X'X
  n <- p - r # Número de filas/columnas dependientes
  m <- n + 1
  A_reducida <- A[-(1:n), -(1:n)] # Extraer submatriz no singular
  D_inv <- solve(A_reducida) # Inversa de la parte no singular
  G <- matrix(0, nrow = p, ncol = p) # Inicializar matriz G con ceros
  G[m:p, m:p] <- D_inv # Insertar la inversa en el bloque inferior derecho
  return(G)
}

A <- t(X3) %*% X3; A

library(MASS)
G3 <- invg(A); fractions(G3)

b0_3 <- G3 %*% t(X3) %*% y3; b0_3

# c. ANOVA

n3 <- length(y3); n3

SCReg_3 <- t(y3) %*% X3 %*% G3 %*% t(X3) %*% y3; SCReg_3

SCTotal_3 <- t(y3) %*% y3; SCTotal_3

SCMedia_3 <- n3*mean(y3)^2; SCMedia_3

SCResid_3 <- SCTotal_3 - SCReg_3; SCResid_3

SCReg_m3 <- SCReg_3 - SCMedia_3; SCReg_m3

s2_3 <- SCResid_3/(n3 - qr(X3)$rank); s2_3

CMReg_m3 <- SCReg_m3/(qr(X3)$rank - 1); CMReg_m3

# 3.a.

lambda3_0 <- matrix(c(1, 0.10, 0.20, 0.40, 0.30), nrow = 5); lambda3_0

mediaestimada3 <- t(lambda3_0) %*% b0_3; mediaestimada3

varestimada3 <- s2_3 %*% t(lambda3_0) %*% G3 %*% lambda3_0; varestimada3

# b. La diferencia en el puntaje entre los profesionales y el promedio
#  de los otros tres grupos?

# \alpha_3 = (\alpha_1 + \alpha_2 + \alpha_4) / 3
# 3 alpha3 - alpha1 - alpha2 - alpha4 = 0

lambda3_b <- matrix(c(0, -1/3, -1/3, 1, -1/3), nrow = 5); lambda3_b

# Ver si es estimable

H3 <- G3 %*% t(X3) %*% X3; H3

all.equal(t(lambda3_b) %*% H3, t(lambda3_b))

## Estimación Puntual

estpuntual3_b <- t(lambda3_b) %*% b0_3; estpuntual3_b

var3_b <- s2_3 %*% t(lambda3_b) %*% G3 %*% lambda3_b; var3_b

# c. La diferencia entre un trabajador independiente y un profesional?

lambda3_c <- matrix(c(0, 0, 0, -1, 1), nrow = 5); lambda3_c

# Ver si es estimable

all.equal(t(lambda3_c) %*% H3, t(lambda3_c))

## Estimación Puntual

estpuntual3_c <- t(lambda3_c) %*% b0_3; estpuntual3_c

var3_c <- s2_3 %*% t(lambda3_c) %*% G3 %*% lambda3_c; var3_c

############################################
############################################
# 4.

## \mu_1 = \mu_2 -> \mu_1 - \mu_2 = 0
## \mu_1 = (\mu_3 + \mu_4)/2 -> 2\mu_1 - \mu_3 - \mu_4 = 0
## \mu_1 = (0.1\mu_1 + 0.4\mu_3)/(0.1 + 0.4) -> \mu_1 - \mu_3 = 0

Lambda4 <- matrix(c(0, 1, -1, 0, 0,
                    0, 2, 0, -1, -1,
                    0, 1, 0, -1, 0), nrow = 3, byrow = T); Lambda4

SCH04 <- t(Lambda4 %*% b0_3) %*% solve(Lambda4 %*% G3 %*% t(Lambda4)) %*%
  (Lambda4 %*% b0_3); SCH04

Fcal4 <- SCH04/(3 * s2_3); Fcal4
Fcrit4 <- qf(0.95, 3, n3 - qr(X3)$rank); Fcrit4

############################################
############################################
# 5.




