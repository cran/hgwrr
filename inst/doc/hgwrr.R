## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(hgwrr)

## -----------------------------------------------------------------------------
data(wuhan.hp)
m_sf <- hgwr(
  formula = Price ~ L(d.Water + d.Commercial) + BuildingArea + (Floor.High | group),
  data = wuhan.hp,
  bw = 299
)

## -----------------------------------------------------------------------------
data(multisampling)
m_df <- hgwr(
  formula = y ~ L(g1 + g2) + x1 + (z1 | group),
  data = multisampling$data,
  coords = multisampling$coords
)

## -----------------------------------------------------------------------------
m_df

## -----------------------------------------------------------------------------
summary(m_df)

## -----------------------------------------------------------------------------
summary(m_df, test_hetero = T)

## -----------------------------------------------------------------------------
head(coef(m_df))
head(fitted(m_df))
head(residuals(m_df))

