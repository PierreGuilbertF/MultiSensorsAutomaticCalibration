# MultiSensorsAutomaticCalibration
Automatic calibration of a multi-sensors system

This program propose an automatic multisensor system calibration

Let's suppose that you have a set of sensors (GPS / IMU, lidar, camera, ...)
on a solid frame. The system sensors + solid frames (Solid) constitute a solid
which means that for any couple of points X1, X2 belong to Solid dX1X2 / dt = 0
at any time.

hence, all sensors si have a constant position Tsi and a constant orientation
Rsi in the solid frame coordinate system. We also suppose that you have the trajectory
and the orientation over the time of all sensors relative to their own orientation / position
coordinate system at time t0.

This tool will provide you an M-estimation of the calibration of all sensors using
their trajectory / orientation over the time. To do that, we use the solid constraint
of all sensors trajectory and orientation over the time. these constraints will lead
to a non linear least square cost function solved using a Levenberg-Marquardt algorithm

In the end, the relative Orientation Rsi/sj and position Tsi/sj of the sensors
will be estimated. The algorithm is proof to gaussian noise since it is a
maximum likelihood estimator

For details about the maths, please read equations.pdf in the doc section
