################################################################################
#               INSTALAÇÃO E CARREGAMENTO DE PACOTES NECESSÁRIOS               #
################################################################################

# Pacotes utilizados:
pacotes <- c("plotly","tidyverse","reshape2","knitr","kableExtra",
             "nlme","lmtest","fastDummies","msm","lmeInfo","jtools")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}

# Algoritmo para determinação dos erros-padrão das variâncias no componente de
# efeitos aleatórios (Fonte: MBA Data Science and Analytics - USP/Esalq)

# ATENÇÃO: A função abaixo é plenamente funcional para modelos do tipo HLM2
# e HLM3, desde que estimados pelo pacote nlme

stderr_nlme <- function(model){
  if(base::class(model) != "lme"){
    base::message("Use a lme object model from nlme package")
    stop()}
  resume <- base::summary(model)
  if(base::length(base::names(model$groups))==1){
    m.type <- "HLM2"
  } else if(base::length(base::names(model$groups))==2){
    m.type <- "HLM3"
  }
  if(m.type == "HLM2"){
    vcov_matrix <- model$apVar
    logs_sd_re <- base::attr(vcov_matrix,"Pars")
    if(base::length(logs_sd_re)==2){
      stderr_tau00 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(`RE Components`=base::c("Var(v0j)","Var(e)"),
                                  `Variance Estimatives`= base::c(base::exp(logs_sd_re)[[1]]^2,
                                                                  base::exp(logs_sd_re[[2]])^2),
                                  `Std Err.`=base::c(stderr_tau00,
                                                     stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                            base::exp(logs_sd_re[[2]])^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                                                               base::exp(logs_sd_re[[2]])^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
    else{
      stderr_tau00 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau01 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x4)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(Components=base::c("Var(v0j)","Var(v1j)","Var(e)"),
                                  Estimatives= base::c(base::exp(logs_sd_re)[[1]]^2,
                                                       base::exp(logs_sd_re[[2]])^2,
                                                       base::exp(logs_sd_re[[4]])^2),
                                  Std_Err=base::c(stderr_tau00,
                                                  stderr_tau01,
                                                  stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                            base::exp(logs_sd_re[[2]])^2/stderr_tau01,
                                            base::exp(logs_sd_re[[4]])^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                                                               base::exp(logs_sd_re[[2]])^2/stderr_tau01,
                                                                               base::exp(logs_sd_re[[4]])^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
  }
  if(m.type == "HLM3"){
    vcov_matrix <- model$apVar
    logs_sd_re <-  base::attr(vcov_matrix,"Pars")
    if(base::length(logs_sd_re) == 3){
      stderr_tau_r000 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u000 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x3)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(Components=base::c("Var(t00k)","Var(v0jk)","Var(e)"),
                                  Estimatives=base::c(base::exp(logs_sd_re)[[2]]^2,
                                                      base::exp(logs_sd_re)[[1]]^2,
                                                      base::exp(logs_sd_re)[[3]]^2),
                                  Std_Err=base::c(stderr_tau_u000,
                                                  stderr_tau_r000,
                                                  stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[2]]^2/stderr_tau_u000,
                                            base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                            base::exp(logs_sd_re)[[3]]^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[2]]^2/stderr_tau_u000,
                                                                               base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                                                               base::exp(logs_sd_re)[[3]]^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    } 
    else{
      stderr_tau_r000 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau_r100 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u000 <- msm::deltamethod(~exp(x4)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u100 <- msm::deltamethod(~exp(x5)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x7)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(`RE_Components`=base::c("Var(t00k)","Var(t10k)",
                                                          "Var(v0jk)","Var(v1jk)",
                                                          "Var(e)"),
                                  `Variance Estimatives`=base::c(base::exp(logs_sd_re)[[4]]^2,
                                                                 base::exp(logs_sd_re)[[5]]^2,
                                                                 base::exp(logs_sd_re)[[1]]^2,
                                                                 base::exp(logs_sd_re)[[2]]^2,
                                                                 base::exp(logs_sd_re)[[7]]^2),
                                  `Std Err.`=base::c(stderr_tau_u000,
                                                     stderr_tau_u100,
                                                     stderr_tau_r000,
                                                     stderr_tau_r100,
                                                     stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[4]]^2/stderr_tau_u000,
                                            base::exp(logs_sd_re)[[5]]^2/stderr_tau_u100,
                                            base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                            base::exp(logs_sd_re)[[2]]^2/stderr_tau_r100,
                                            base::exp(logs_sd_re)[[7]]^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[4]]^2/stderr_tau_u000,
                                                                               base::exp(logs_sd_re)[[5]]^2/stderr_tau_u100,
                                                                               base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                                                               base::exp(logs_sd_re)[[2]]^2/stderr_tau_r100,
                                                                               base::exp(logs_sd_re)[[7]]^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
  }
}

# Carregando a base de dados
load('tempo_estudante_escola.RData')

# Visualizando a base de dados
tempo_estudante_escola %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                font_size = 12, full_width = TRUE)

glimpse(tempo_estudante_escola)

# Analisando estatísticas descritivas e univariadas:
summary(tempo_estudante_escola)

# Balanceamento dos dados em relação à quantidade de alunos por período 
# analisado (mês)
tempo_estudante_escola %>% 
  rename(Mês = 3, 
         `Quantidade de alunos` = 2) %>% 
  group_by(Mês) %>% 
  summarise(`Quantidade de Alunos` = n()) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                full_width = TRUE, font_size = 12)

# Balanceamento dos dados em relação à quantidade de alunos por escola
tempo_estudante_escola %>% 
  rename(Escola = 1, 
         `Quantidade de Alunos` = 2) %>% 
  group_by(Escola) %>% 
  summarise(`Quantidade de Alunos` = n()/4) %>%
  kable() %>%                                   
  kable_styling(bootstrap_options = 'striped', 
                full_width = TRUE, font_size = 12)

# Note que só podemos usar o artifício de dividir o número de estudantes por 4
# porque sabemos que a quantidade de alunos é igual nos 4 períodos de tempo 
# analisados (amostragem balanceada - alunos x tempo)

# Evolução temporal média do desempenho escolar dos estudantes (utilizando 
# ajuste linear apenas)

ggplotly(
  tempo_estudante_escola %>% 
    ggplot(aes(x = mes, y = desempenho, label = estudante))+
    geom_point(color = 'aquamarine4', alpha = .1) +
    geom_smooth(color = 'darkorchid4', method = 'lm', se = F, size = 2) +
    labs(x = 'Mês', y = "Desempenho escolar") +
    theme_bw()
)

# Kerner Density Estimation (KDE) - plotando a função densidade de probabilidade
# da variável dependente (desempenho) com histograma

ggplotly(
  tempo_estudante_escola %>% 
    ggplot() +
    geom_histogram(aes(x = desempenho, y = ..density..), color = 'black', 
                   fill = 'aquamarine4')+
    geom_density(aes(x = desempenho), color = 'black', size = 1) + 
    labs(x = "Desempenho", y = "Densidade") +
    theme_bw()
)

# Kerner Density Estimation (KDE) - plotando a função densidade de probabilidade
# da variável dependente (desempenho) por escola (incluindo a variável escola
# dentro do aesthetic do geom_density)

ggplotly(
  tempo_estudante_escola %>% 
    ggplot() +
    geom_histogram(aes(x = desempenho, y = ..density..), color = 'black', 
                   fill = 'aquamarine4', alpha = .5)+
    geom_density(aes(x = desempenho, color = escola, fill = escola), 
                 alpha = .2) + 
    labs(x = "Desempenho", y = "Densidade") +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    theme_bw()
)


# Kerner Density Estimation (KDE) - plotando a função densidade de probabilidade
# da variável dependente (desempenho) com histograma e por escola (função 
# facet_wrap)

tempo_estudante_escola %>% 
  group_by(escola) %>% 
  mutate(linhas = 1:n()) %>% 
  mutate(x = unlist(density(desempenho, n = max(linhas))['x']),
         y = unlist(density(desempenho, n = max(linhas))['y'])) %>% 
  ggplot() +
  geom_area(aes(x = x, y = y, group = escola, fill = escola), color = 'black', 
            alpha = .3) + 
  geom_histogram(aes(x = desempenho, y = ..density.., fill = escola), 
                 color = 'black', alpha = .1, position = 'identity') +
  facet_wrap(~escola) + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  theme_bw()


  
