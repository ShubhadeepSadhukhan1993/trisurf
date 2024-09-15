/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _TIMESTEP_H
#define _TIMESTEP_H
ts_bool single_timestep(ts_vesicle *vesicle, ts_double *vmsr, ts_double *bfsr);
ts_bool run_simulation(ts_vesicle *vesicle, ts_uint mcsweeps, ts_uint inititer, ts_uint iterations, ts_uint start_simulation);

void calculate_proj_max_min(ts_vesicle *vesicle);
void UCSP_step(ts_vesicle *vesicle);
double error_check(ts_vesicle *vesicle);
void calculate_vesicle_dir(ts_vesicle *vesicle, double dir[3]);

void calculate_velocity_mag(ts_vesicle *vesicle, double *vmag);
void calculate_xback(ts_vesicle *vesicle);
#endif
