# Copyright (c) 2025 Bart van de Lint
# SPDX-License-Identifier: MPL-2.0

"""
Real-time 3D visualization with keyboard steering.

Controls:
  Arrow keys to steer, ESC to stop.
"""

using GLMakie
using SymbolicAWEModels, VortexStepMethod, KiteUtils
using SymbolicAWEModels: init!, next_step!, update_sys_state!, log!, calc_steady_torque, save_log, Logger
using Printf

# Parameters
dt = 0.05
total_time = 60.0
vsm_interval = 3
realtime_factor = 1.0
plot_interval = 1
vector_scale = 1.0
steering_magnitude = 5.0
record_video = false
output_filename = "data/kite_simulation.mp4"
framerate = 20

# Initialize model
set_data_path("data/2plate_kite")
struc_yaml = joinpath(get_data_path(),
                      "quat_struc_geometry.yaml")
aero_yaml = joinpath(get_data_path(), "aero_geometry.yaml")
update_aero_yaml_from_struc_yaml!(struc_yaml, aero_yaml)

set = Settings("system.yaml")
set.profile_law = 3
vsm_set = VortexStepMethod.VSMSettings(
    joinpath(get_data_path(), "vsm_settings.yaml");
    data_prefix=false)
sys = load_sys_struct_from_yaml(struc_yaml;
    system_name="2plate_kite", set, vsm_set)
sam = SymbolicAWEModel(set, sys)
init!(sam; remake=false)
find_steady_state!(sam)

sys_struct = sam.sys_struct
steps = Int(round(total_time / dt))

# Create plot
scene = plot(sys_struct; vector_scale, size=(1400, 900))
display(scene)

progress_text = Observable("Progress: 0%")
text!(scene, progress_text,
      position=Point2f(1380, 40), space=:pixel,
      fontsize=20, color=:black, align=(:right, :top))
text!(scene, "Arrow keys to steer, ESC to stop",
      position=Point2f(20, 40), space=:pixel,
      fontsize=16, color=:darkblue, align=(:left, :top))

# Keyboard control
current_steering = Observable([0.0, 0.0, 0.0])
stop_simulation = Ref(false)

on(events(scene).keyboardbutton) do event
    mag = steering_magnitude
    if event.action in (Keyboard.press, Keyboard.repeat)
        if event.key == Keyboard.left
            current_steering[] = [0.0, -mag, mag]
        elseif event.key == Keyboard.right
            current_steering[] = [0.0, mag, -mag]
        elseif event.key == Keyboard.down
            current_steering[] = [0.0, -mag, -mag]
        elseif event.key == Keyboard.up
            current_steering[] = [0.0, mag, mag]
        elseif event.key == Keyboard.escape
            stop_simulation[] = true
        end
    elseif event.action == Keyboard.release &&
           event.key != Keyboard.escape
        current_steering[] = [0.0, 0.0, 0.0]
    end
end

function run_realtime!(sam, sys_struct, scene)
    logger = Logger(sam, steps + 1)
    sys_state = SysState(sam)
    sys_state.time = 0.0
    log!(logger, sys_state)

    steady_torque =
        calc_steady_torque(sam)
    torque_damp = 0.9
    start_time = time()
    sim_time = 0.0
    io = record_video ?
        VideoStream(scene; framerate) : nothing

    for step in 1:steps
        stop_simulation[] && break
        t = step * dt
        target_elapsed = t / realtime_factor

        steady_torque = torque_damp * steady_torque +
            (1 - torque_damp) *
            calc_steady_torque(sam)
        control = steady_torque .+ current_steering[]

        t0 = time()
        next_step!(sam; set_values=control,
                   dt, vsm_interval)
        sim_time += time() - t0

        update_sys_state!(
            sys_state, sam)
        sys_state.time = t
        log!(logger, sys_state)

        if step % plot_interval == 0
            plot!(sys_struct; vector_scale)
            progress_text[] = @sprintf(
                "Progress: %d%%",
                round(Int, 100 * step / steps))
            record_video && recordframe!(io)
            sleep(0.001)
        end

        actual_elapsed = time() - start_time
        sleep(max(0.0, target_elapsed - actual_elapsed))

        if step % (steps ÷ 10) == 0
            @printf("  %.0f%% (t=%.1fs)\n",
                    100 * step / steps, t)
        end
    end

    record_video && save(output_filename, io)

    elapsed = time() - start_time
    runtime = round(elapsed; digits=2)
    sim_speedup = round(total_time / sim_time; digits=2)
    @info "Done. Runtime: " runtime "Simulation speedup: " sim_speedup

    save_log(logger,
                               "tmp_realtime_run")
    lg = load_log("tmp_realtime_run")
    display(plot(sam.sys_struct, lg))
end

run_realtime!(sam, sys_struct, scene)
