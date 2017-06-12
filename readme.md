# Multi-Trajectory Pose Correspondence in Pose-Graphs

This code considers the problem of finding pose matches between trajectories
of multiple robots in their respective coordinate frames or equivalent matches
between trajectories obtained during different sessions. Pose correspondences
between trajectories are mediated by common landmarks represented in a
topological map lacking distinct metric coordinates. Despite such lack of
explicit metric level associations, we mine preliminary pose level
correspondences between trajectories through a novel multi-scale heat-kernel 
descriptor and correspondence graph framework. These serve as an improved 
initialization for ICP (Iterative Closest Point) to yield dense pose 
correspondences.
 


## Getting Started

Run in matlab

```
mainRUN
```

## Authors

* **Sayantan Datta** - *Initial work* - [kenzo450D](https://github.com/kenzo450D)


