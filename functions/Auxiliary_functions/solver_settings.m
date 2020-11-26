function solver_settings =solver_settings() 
solver_settings.iterations_max          = 400                       ;
common_start                            = 10;
solver_settings.epsilon                 = 1e-3                      ;
solver_settings.epsilon_feas            = 1e-2                      ;
solver_settings.t_stop                  = 1000                     ;
solver_settings.records_num             = 100                       ;
solver_settings.nu                      = 10                       ; 
solver_settings.alpha                   = 0.99                      ;
solver_settings.beta                    = 0.05                      ;
solver_settings.theta                   = 0.95                      ;
solver_settings.s_efficient             = 1e-3                      ;
solver_settings.lambda_correction       = 0.5                       ;
solver_settings.t_0                     = 1                         ;
solver_settings.x_0                     = 0                        ;
solver_settings.lambda_0                = common_start              ; 
solver_settings.gamma_0                 = common_start              ; 
solver_settings.mu_0                    = common_start              ; 
solver_settings.slack_0                 = .99                       ;


end

