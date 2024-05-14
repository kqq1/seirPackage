#' SEI1I2R传染病模型
#'
#' 该函数实现了一个包含易感者(S)、暴露者(E)、首次感染者(I1)、二次感染者(I2)和康复者(R)的传染病模型。
#' @param t 时间变量。
#' @param y 状态变量的向量，包含S, E, I1, I2, R。
#' @param params 参数列表，包含总人口数(N)、感染率(beta)、从易感者到暴露者的速率(sigma)、首次感染到康复的速率(gamma1)、二次感染到康复的速率(gamma2)和自然死亡率(mu)。
#' @return 一个包含状态变量导数的向量。
#' @export
# 定义模型参数
params <- list(N = 1000, beta = 0.5, sigma = 0.1, gamma1 = 0.2, gamma2 = 0.3, mu = 0.01)

# 初始条件
yinit <- c(S = 990, E = 10, I1 = 0, I2 = 0, R = 0)

# 时间序列
times <- seq(0, 200, by = 0.1)

# SEI1I2R传染病模型

# 定义模型参数
params <- list(N = 1000, beta = 0.5, sigma = 0.1, gamma1 = 0.2, gamma2 = 0.3, mu = 0.01)

# 初始条件
yinit <- c(S = 990, E = 10, I1 = 0, I2 = 0, R = 0)

# 时间序列
times <- seq(0, 200, by = 0.1)

# 定义模型函数
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # 计算导数
    dS <- -beta * S * (I1 + I2) / N
    dE <- beta * S * (I1 + I2) / N - sigma * E
    dI1 <- sigma * E - gamma1 * I1 - mu * I1
    dI2 <- gamma1 * I1 - gamma2 * I2 - mu * I2
    dR <- gamma2 * I2

    # 返回导数列表
    list(c(dS, dE, dI1, dI2, dR))
  })
}

# 运行模型模拟
out <- ode(y = yinit, times = times, func = model, parms = params)

# 绘制SEIR图
seirplot <- function(times, out, N) {
  library(ggplot2)

  # 将out转换为数据框
  out_df <- as.data.frame(out)

  # 计算比率
  out_df$S_ratio <- out_df$S / N
  out_df$E_ratio <- out_df$E / N
  out_df$I1_ratio <- out_df$I1 / N
  out_df$I2_ratio <- out_df$I2 / N
  out_df$R_ratio <- out_df$R / N

  plot <- ggplot(data = out_df, aes(x = time)) +
    geom_line(aes(y = S_ratio, color = "S")) +
    geom_line(aes(y = E_ratio, color = "E")) +
    geom_line(aes(y = I1_ratio, color = "I1")) +
    geom_line(aes(y = I2_ratio, color = "I2")) +
    geom_line(aes(y = R_ratio, color = "R")) +
    scale_y_continuous(name = "Ratio", limits = c(0, 1), breaks = seq(0, 1, by = 0.24)) +
    labs(title = "SEI1I2R Model", x = "Time", y = "Ratio", color = "Legend") +
    scale_color_manual(values = c("S" = "blue", "E" = "red", "I1" = "green", "I2" = "orange", "R" = "purple"))

  print(plot)
}

# 定义总体人数N
N <- 1000

# 重新运行模拟
out <- ode(y = yinit, times = times, func = model, parms = params)
# 绘制SEIR图
seirplot(times, out, N)
