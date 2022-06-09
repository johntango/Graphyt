 

#ifndef BOUNDARYCOND_H
#define BOUNDARYCOND_H
// Enum for BC types
typedef struct{
    enum function{
        none,
        grid_vx,
        grid_vy,
        grid_vz,
        grid_temp,
        grid_tflux_x,
        grid_tflux_y,
        grid_tflux_z,
        grid_fx,
        grid_fy,
        grid_fz,
        grid_vx_vy,
        grid_vx_vz,
        grid_vy_vz,
        grid_fx_fy,
        grid_fy_fz,
        grid_fx_fz,
        grid_fx_fy_fz,
        grid_vx_vy_vz,
        grid_vx_temp,
        grid_vy_temp,
        grid_vz_temp,
        grid_fx_temp,
        grid_fy_temp,
        grid_fz_temp,
        grid_vx_vy_temp,
        grid_vx_vz_temp,
        grid_vy_vz_temp,
        grid_fx_fy_temp,
        grid_fy_fz_temp,
        grid_fx_fz_temp,
        grid_fx_fy_fz_temp,
        grid_vx_vy_vz_temp,
        mp_vx,
        mp_vy,
        mp_vz,
        mp_temp,
        mp_fx,
        mp_fy,
        mp_fz,
        mp_vx_vy,
        mp_vx_vz,
        mp_vy_vz,
        mp_fx_fy,
        mp_fy_fz,
        mp_fx_fz,
        mp_fx_fy_fz,
        mp_vx_vy_vz,
        mp_vx_temp,
        mp_vy_temp,
        mp_vz_temp,
        mp_fx_temp,
        mp_fy_temp,
        mp_fz_temp,
        mp_vx_vy_temp,
        mp_vx_vz_temp,
        mp_vy_vz_temp,
        mp_fx_fy_temp,
        mp_fy_fz_temp,
        mp_fx_fz_temp,
        mp_fx_fy_fz_temp,
        mp_vx_vy_vz_temp,
        mp_stress,
        mp_stress_xx_yy,
        mpstrain,
        mp_vxsin,
        mp_vysin,
        mp_syysin
    };
} BCTypes;





#endif // BOUNDARYCOND_H
