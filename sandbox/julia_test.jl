
using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")

using CSV, DataFrames

pwd()

seed_dat = DataFrame(CSV.File("./data/seed_effects_data_1.csv"))


names(seed_dat)

first(seed_dat, 10)

last(seed_dat, 10)

## random forest 

Pkg.add("ScikitLearn")
Pkg.add("BenchmarkTools")

using CSV
using DataFrames
using ScikitLearn: @sk_import, fit!, predict
@sk_import ensemble: RandomForestRegressor
using ScikitLearn.GridSearch: RandomizedSearchCV
using BenchmarkTools

seed_dat = CSV.read("./data/seed_effects_data_1.csv", DataFrames.DataFrame)

outcome = seed_dat[:,"pree_acog"]
exposure = seed_dat[:,"fv_totdens_2_5"]
covariates = select!(seed_dat, Not([:"pree_acog", :"fv_totdens_2_5"]))

mod = RandomForestRegressor(max_leaf_nodes=2)
param_dist = Dict("n_estimators"=>[50 , 100, 200, 300],
                  "max_depth"=> [3, 5, 6 ,8 , 9 ,10])

model = RandomizedSearchCV(mod, param_dist, n_iter=10, cv=5, n_jobs=1)

@btime fit!(model, Matrix(DataFrames(covariates)), Matrix(DataFrames(outcome)))