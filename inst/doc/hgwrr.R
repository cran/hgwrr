## ---- include = FALSE---------------------------------------------------------
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
m_sf

## -----------------------------------------------------------------------------
summary(m_sf)

## -----------------------------------------------------------------------------
head(coef(m_sf))
head(fitted(m_sf))
head(residuals(m_sf))

