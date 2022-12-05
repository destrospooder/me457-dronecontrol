# ME457 - Drone Control (Fall 2022)

Repository containing programming assignments for ME457: Drone Control, taken during Fall 2022. Project collaborators are Sohaib Bhatti, Kyle Deolall, and Ariel Tamayev. Probably not a good idea to keep this public during the semester because of mucky no-good code stealers, but what can you do? Oh well.

**Instructor:** Dirk Luchtenburg<br/>
**Email:** dirk.luchtenburg@cooper.edu<br/>

**Textbook:** "[Small Unmanned Aircraft: Theory and Practice](https://github.com/randybeard/uavbook)" (2) by Randal Beard and Timothy W. McLain

**Description:**<br/>
This course prepares students to do research in the rapidly evolving of field of autonomous navigation, guidance, and control of unmanned air vehicles (UAVs). In particular, students will learn about key concepts from rigid-body dynamics, aerodynamics, feedback control, and state estimation using sensors, to maneuver through obstacles. Traditional homework assignments are replaced with a semester-long simulation software development project in Python. Techniques developed will be applied in the form of student design projects.

---

Log:

*10/9 - Chapter 2, Chapter 3 sort of complete, weird bug where plane keeps going north* -BA

*10/10 - All good! (Findings from Presentation 1 - Look into gyroscopic coupling. The reason it's rolling eccentrically is because of the J_xz; if you set it equal to 0, it ends up rolling symmetrically along the x-axis. MAKE GRAPHS BIGGER!)* -BA
 
*10/23 - All not good (mostly good). Need to finish `mav_dynamics`, there's some wind stuff there I don't feel like finishing (sick of looking at code), and some thrust stuff too. Wind simulation needs transfer function stuff as well. Updated Aerosonde with E2 Aerosonde Stuff (p. 276 of Beard). The wind simulation parameters Va...Lw, and all the sigmas are from the book, not sure if they're correct (4.4). Also Sohaib says there's a fundamental dynamics problem from week 1 (10/10) which still needs to be solved I believe.* -AT

*10/24 - Got the sim working so that's something. The dynamics are still screwy but we'll get people to look through it.* -BA

*11/16 - Unit tests all passed! PID control for DC motor included.* -SB

*12/5 - Call man fill tuhr.* -BA
