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

# Explorando visualmente a evolução temporal dos 50 primeiros estudantes da BD (
# limitação do visual do gráfico)
glimpse(tempo_estudante_escola)

tempo_estudante_escola %>%
  filter(estudante %in% 1:50) %>% 
  ggplot(aes(x = factor(mes), y = desempenho, group = estudante)) +
  geom_point(aes(color = estudante), size = 2) +
  geom_line(aes(color = estudante), size = 1) +
  guides(color = "none") +
  labs(x = "Mês", y = "Desempenho Escolar") +
  scale_color_viridis_d() +
  theme_bw()
  

# Gráfico da evolução temporal (medidas repetidas / longitudinal) do desempenho
# médio dos alunos por escola

ggplotly(
  tempo_estudante_escola %>% 
    ggplot(aes(x = mes, y = desempenho, color = escola)) +
    geom_point(size = 2.5, alpha = .3) + 
    geom_smooth(aes(group = escola), method = 'lm', se = F, formula = y ~ x) +
    labs(x = "Mês", y = "Desempenho Escolar", color = "Escolas:") +
    guides(color = 'none') +
    scale_color_viridis_d() +
    theme_bw()
)


# Estimação do modelo HLM3 nulo - 1º passo do step up -------------------------

modelo_nulo_hlm3 <- lme(fixed = desempenho ~ 1, 
                        random = list(escola = ~ 1, estudante = ~ 1), 
                        data = tempo_estudante_escola, 
                        method = 'REML')

# Parâmetros do modelo
summary(modelo_nulo_hlm3)

# Cálculo dos erros padrão do modelo por meio da função stderr_nlme (Fonte no 
# início do script) 
stderr_nlme(modelo_nulo_hlm3)


# Comparando o modelo multinível nulo com OLS nulo ----------------------------

modelo_ols_nulo <- lm(formula = desempenho ~ 1, 
                      data = tempo_estudante_escola)

# Parâmetros do modelo
summary(modelo_ols_nulo)

# Comparando os LogLiks (função do pacote lmtest)
lrtest(modelo_nulo_hlm3, modelo_ols_nulo)

# Gráfico apresentando a evolução do step up
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM3_Nulo = logLik(modelo_nulo_hlm3)) %>%
  rename(`OLS Nulo` = 1,
         `HLM3 Nulo` = 2) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.1, color = "white", 
             size = 7) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("aquamarine4","darkorchid4")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())


# Estimação do modelo HLM3 com tendênca linear e Interceptos e Inclinações
# aleatórias ------------------------------------------------------------------

modelo_intercept_inclin_hlm3 <- lme(fixed = desempenho ~ mes, 
                                    random = list(escola = ~ mes,
                                                  estudante = ~ mes),
                                    data = tempo_estudante_escola, 
                                    method = "REML")

# Parâmetros do modelo
summary(modelo_intercept_inclin_hlm3)

# Cálculo dos erros padrão do modelo por meio da função stderr_nlme (Fonte no 
# início do script) 
stderr_nlme(modelo_intercept_inclin_hlm3)

# Comparando os LogLiks (função do pacote lmtest)
lrtest(modelo_intercept_inclin_hlm3, modelo_nulo_hlm3)

# Gráfico apresentando a evolução do step up
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM3_Nulo = logLik(modelo_nulo_hlm3),
           HLM3_Intercept_Incli = logLik(modelo_intercept_inclin_hlm3)) %>%
  rename(`OLS Nulo` = 1,
         `HLM3 Nulo` = 2,
         `HLM3 c/ Intercept. e Inclin. aleatórios` = 3) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.1, color = "white", 
             size = 5) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("aquamarine4","darkorchid4", 'darkred')) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())


# Estimação do modelo HLM3 com tendência linear, interceptos e inclinações 
# aleatórias e variáveis de nível 2 (ativ - dummie relativa a realização ou não
# de atividades extras) e nível 3 (texp - tempo de experiência médio dos 
# professores de determinada escola) ------------------------------------------

modelo_completo_hlm3 <- lme(fixed = desempenho ~ mes + ativ + texp + 
                              ativ:mes + texp:mes, 
                            random = list(escola = ~ mes, 
                                          estudante = ~ mes), 
                            data = tempo_estudante_escola, 
                            method = 'REML')

# Parâmetros do modelo
summary(modelo_completo_hlm3)

# Cálculo dos erros padrão do modelo por meio da função stderr_nlme (Fonte no 
# início do script) 
stderr_nlme(modelo_completo_hlm3)

# Comparando os LogLiks (função do pacote lmtest)
lrtest(modelo_completo_hlm3,modelo_intercept_inclin_hlm3)

# Gráfico apresentando a evolução do step up
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM3_Nulo = logLik(modelo_nulo_hlm3),
           HLM3_Intercept_Incli = logLik(modelo_intercept_inclin_hlm3),
           HLM3_completo = logLik(modelo_completo_hlm3)) %>%
  rename(`OLS Nulo` = 1,
         `HLM3 Nulo` = 2,
         `HLM3 c/ Intercept. e Inclin. aleatórios` = 3,
         `HLM3 completo` = 4) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.1, color = "white", 
             size = 5) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("aquamarine4","darkorchid4", 'darkred',
                               'darkblue')) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

# Apresentando os valores dos efeitos aleatórios de intercepto e inclinação do
# nível estudante e escola, respectivamente: v0jk, v1jk, t00k, t10k 

# Nível Estudante
random.effects(modelo_completo_hlm3)[['estudante']] %>% 
  rename(v0jk = 1, v1jk = 2) %>% 
  rownames_to_column("Estudante") %>% 
  mutate(Estudante = gsub("^.*?\\/", "", Estudante)) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                full_width = T, font_size = 12)

# Nível Escola
random.effects(modelo_completo_hlm3)[['escola']] %>% 
  rename(t00k = 1, t10k = 2) %>% 
  rownames_to_column("Escola") %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                full_width = T, font_size = 12)

# Plotando o comportamento de v0jk (interceptos aleatórios no nível estudante)
# para fins didáticos (gráfico poluído devido a excesso de dados)

ggplotly(
  random.effects(modelo_completo_hlm3)[['estudante']] %>% 
    rename(v0jk = 1, v1jk = 2) %>% 
    rownames_to_column("Estudante") %>% 
    mutate(Estudante = gsub("^.*?\\/", "", Estudante)) %>% 
    ggplot(aes(x = fct_rev(Estudante), y = v0jk, label = Estudante)) +
    geom_bar(stat = 'identity', color = 'aquamarine4') + 
    coord_flip() +
    labs(x = "Estudante", y = "v0jk") +
    theme(legend.title = element_blank(),
          panel.background = element_rect('white'),
          legend.position = 'none', 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
)

# Plotando o comportamento de v1jk (inclinações aleatórias no nível estudante)
# para fins didáticos (gráfico poluído devido a excesso de dados)

ggplotly(
  random.effects(modelo_completo_hlm3)[['estudante']] %>% 
    rename(v0jk = 1, v1jk = 2) %>% 
    rownames_to_column("Estudante") %>% 
    mutate(Estudante = gsub("^.*?\\/", "", Estudante)) %>% 
    ggplot(aes(x = fct_rev(Estudante), y = v1jk, label = Estudante)) +
    geom_bar(stat = 'identity', color = 'darkorchid4') + 
    coord_flip() +
    labs(x = "Estudante", y = "v1jk") +
    theme(legend.title = element_blank(),
          panel.background = element_rect('white'),
          legend.position = 'none', 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
)

# Plotando o comportamento de t00k (interceptos aleatórios no nível escola)
# para fins didáticos

random.effects(modelo_completo_hlm3)[["escola"]] %>% 
  rename(t00k = 1,
         t10k = 2) %>% 
  rownames_to_column("Escola") %>% 
  mutate(color_t00k = ifelse(t00k < 0, "A", "B"),
         color_t10k = ifelse(t10k < 0, "A", "B"),
         hjust_t00k = ifelse(t00k > 0, 1.15, -0.15),
         hjust_t10k = ifelse(t10k > 0, 1.15, -0.15)) %>% 
  arrange(Escola) %>%
  ggplot(aes(label = round(t00k, digits = 3), 
             hjust = hjust_t00k)) +
  geom_bar(aes(x = fct_rev(Escola), y = t00k, fill = color_t00k),
           stat = "identity", color = "black") +
  geom_text(aes(x = Escola, y = 0), size = 4.1, color = "black") +
  coord_flip() +
  labs(x = "Escola",
       y = expression(t[0][0][k])) +
  scale_fill_manual(values = c("darkred","darkgreen")) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey95"),
        legend.position = "none")

# Plotando o comportamento de t10k (inclinações aleatórias no nível escola)
# para fins didáticos

random.effects(modelo_completo_hlm3)[["escola"]] %>% 
  rename(t00k = 1,
         t10k = 2) %>% 
  rownames_to_column("Escola") %>% 
  mutate(color_t00k = ifelse(t00k < 0, "A", "B"),
         color_t10k = ifelse(t10k < 0, "A", "B"),
         hjust_t00k = ifelse(t00k > 0, 1.15, -0.15),
         hjust_t10k = ifelse(t10k > 0, 1.15, -0.15)) %>% 
  arrange(Escola) %>% 
  ggplot(aes(label = round(t10k, digits = 3), 
             hjust = hjust_t10k)) +
  geom_bar(aes(x = fct_rev(Escola), y = t10k, fill = color_t10k),
           stat = "identity", color = "black") +
  geom_text(aes(x = Escola, y = 0), size = 4.1, color = "black") +
  coord_flip() +
  labs(x = "Escola",
       y = expression(t[1][0][k])) +
  scale_fill_manual(values = c("darkred","darkgreen")) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey95"),
        legend.position = "none")

# Fitted values do modelo
predict(modelo_completo_hlm3, level = 0:2) %>% 
  mutate(estudante = gsub('^.*?\\/', '', estudante),
         estudante = as.factor(as.numeric(estudante)),
         escola = as.factor(as.numeric(escola)),
         mes = tempo_estudante_escola$mes,
         desempenho = tempo_estudante_escola$desempenho,
         etjk = resid(modelo_completo_hlm3)) %>%
  rename("fitted_fixed" = 3, 
         "fitted_escola" = 4, 
         "fitted_estudante" = 5) %>% 
  select(escola, estudante, mes, desempenho, everything()) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = 'striped', 
                full_width = T, font_size = 12)

# Exemplos de predições com o modelo estimado ---------------------------------

# Exemplo: Quais os valores previstos de desempenho escolar no primeiro mês 
# (mes = 1) para o estudante "1" da escola "1", sabendo-se que esta escola 
# oferece tempo médio de experiência de seus professores igual a 2 anos?
glimpse(tempo_estudante_escola)

predict(object = modelo_completo_hlm3, level = 0:2,
        newdata = data.frame(mes = 1, estudante = "1", escola = "1",texp = 2,
                             ativ = c("sim", "não")))

# Fitted values (desempenho previsto) dos 47 alunos da primeira escola
predict(modelo_completo_hlm3, level = 0:2) %>% 
  mutate(estudante = gsub("^.*?\\/","",estudante),
         estudante = as.factor(as.numeric(estudante)),
         mes = tempo_estudante_escola$mes) %>% 
  rename(fitted_fixed = 3,
         fitted_escola = 4,
         fitted_estudante = 5) %>% 
  right_join(tempo_estudante_escola, 
             by = c("escola","estudante","mes")) %>%
  filter(estudante %in% 1:47) %>%
  ggplot(aes(x = mes, y = fitted_estudante, color = estudante)) +
  geom_smooth(method = "lm", formula = y ~ x, se = F, size = 1.3) +
  geom_point(size = 4, alpha = 0.4) +
  guides(color = F) +
  scale_colour_viridis_d() +
  labs(x = "Mês",
       y = "Desempenho Escolar") +
  theme_bw()


# Bônus - comparando o modelo multinível com o modelo OLS utilizando dummies
# nas variáveis que fazem parte do contexto multinível (superficialmente, uma
# comparação do que seria considerar esses níveis em uma modelagem tradicional
# do tipo GLM, por OLS)

# Serão utilizadas para o modelo OLS as mesma variáveis utilizadas no multinível

# Aplicando o procedimento n-1 dummies na base de dados
base_dummizada <- dummy_cols(.data = tempo_estudante_escola, 
                             select_columns = 'escola', 
                             remove_most_frequent_dummy = T, 
                             remove_selected_columns = T)

# Estimado a OLS
modelo_ols_dummies <- lm(formula = desempenho ~ . + mes:texp + mes:ativ -
                           estudante, data = base_dummizada)

#Parâmetros do modelo ols_final
summary(modelo_ols_dummies)

# Aplicando o procedimento Step wise para retirar as variáveis não 
# estatisticamente significantes
modelo_ols_dummies_step <- step(object = modelo_ols_dummies, )

# Parâmetros do modelo OLS final - Step Wise
summary(modelo_ols_dummies_step)

# Comparando os LogLiks dos dois modelos
lrtest(modelo_completo_hlm3, modelo_ols_dummies_step)

