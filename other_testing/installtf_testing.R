library(reticulate)
use_condaenv("r-tensorflow")
import("tensorflow", as = "tf")

###
library(reticulate)
use_python('~/.virtualenvs/r-tensorflow/Scripts/python.exe')
library(keras)

model <- keras_model_sequential() |>
  layer_dense(units = 100, activation = 'relu', input_shape = 1000) |>
  layer_dense(units = 32, activation = 'relu') |>
  layer_dense(units = 6, activation = 'linear')

###
Sys.setenv(RETICULATE_PYTHON = '~/.virtualenvs/r-tensorflow/Scripts')
library(tensorflow)

###
library(tensorflow)
use_condaenv("r-tensorflow")
