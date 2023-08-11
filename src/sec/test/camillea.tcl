# test for camellia die
source "helpers.tcl"

read_lef asap7_tech_4x_201209.lef
read_lef asap7sc7p5t_28_SL_4x_220121a.lef
read_lef asap7sc7p5t_28_L_4x_220121a.lef

read_def "camellia.def"

sec_evaluate
