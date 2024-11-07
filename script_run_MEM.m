data_in = 'data/GMQMLADKNU/in/data_tapping_start_average_240724_1624.mat';
options_in = 'wMEM_options.json';


profile on
run_MEM(data_in, 'cMEM_options.json')
                                         
profile viewer
