"""
mavsimPy
    - Chapter 4 assignment for Beard & McLain, PUP, 2012
    - Update history:
        12/27/2018 - RWB
        1/17/2019 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
import parameters.simulation_parameters as SIM

from chap2.mav_viewer import MavViewer
from chap3.data_viewer import DataViewer
from chap4.mav_dynamics import MavDynamics
from chap4.wind_simulation import WindSimulation
from message_types.msg_delta import MsgDelta

# initialize the visualization
VIDEO = False  # True==write video, False==don't write video
mav_view = MavViewer()  # initialize the mav viewer
if VIDEO is True:
    from chap2.video_writer import VideoWriter
    video = VideoWriter(video_name="chap4_video.avi",
                        bounding_box=(0, 0, 1000, 1000),
                        output_rate=SIM.ts_video)

# initialize elements of the architecture
wind = WindSimulation(SIM.ts_simulation)
mav = MavDynamics(SIM.ts_simulation)
delta = MsgDelta()

# initialize the simulation time
sim_time = 0
end_time = 1000
plot_time = sim_time
refr_rate = 0.1
time_step = 0.02

# main simulation loop
print("Press Command-Q to exit...")
while sim_time < end_time:
    # -------set control surfaces-------------
    delta.elevator = 0.5
    delta.aileron = 0.5
    delta.rudder = 0.5
    delta.throttle = 0.5

    # -------physical system-------------
    current_wind = wind.update()  # get the new wind vector
    mav.update(delta, current_wind)  # propagate the MAV dynamics

    # -------update viewer-------------
    if sim_time-plot_time > refr_rate:
        mav_view.update(mav.true_state)  # plot body of MAV
        plot_time = 500
    if VIDEO is True:
        video.update(sim_time)

    # -------increment time-------------
    sim_time += time_step

if VIDEO is True:
    video.close()

