# first, you need to prepare your "grid data", see _utils/read_grid_data for more details
butane_2d_log_prob_example.txt
pentane_3d_log_prob_example.txt

# second, you need to convert your grid data into VTI files
python make_ttk_input.py butane_2d_log_prob_example.txt butane_2d_log_prob_example.vti
python make_ttk_input.py pentane_3d_log_prob_example.txt pentane_3d_log_prob_example.vti

# third, you can run TTK
python run_ttk.py butane_2d_log_prob_example.vti butane_results 1 1 0
python run_ttk.py pentane_3d_log_prob_example.vti pentane_results 1.1 2 3

# finally, you can visualize your results
python plot_ttk_results.py butane_2d_log_prob_example.txt butane_results
python plot_ttk_results.py pentane_3d_log_prob_example.txt pentane_results