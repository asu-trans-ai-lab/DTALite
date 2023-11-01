# DTALite

This repository maintains the source code and the Windows-based release for DTALite+NeXTA
package. DTALite is an open-source AMS library for efficiently macroscopic and
mesoscopic traffic assignment based on [General Modeling Network
Specification](https://github.com/zephyr-data-specs/GMNS) (GMNS) format. NeXTA is developed through the support of a [FHWA
study](https://www.fhwa.dot.gov/publications/research/operations/13036/004.cfm) as a visualization tool for transportation Analysis, Modeling, and Simulation(AMS).

**Some Remarks:**

1. [Path4GMNS](https://github.com/jdlph/Path4GMNS) provides an Python API to the portable DTALite running across Windows,
   Linux, and macOS (including both x86_64 and arm64 architectures).
2. Please download the latest [DTALite-NeXTA package](release/), which follows the calendar versioning scheme.

## Getting Started

- [Self-learning documents](docs/self_learning)
- [User guide for NeXTA+QGIS using GMNS data format](docs/user_guide/1_QGIS_NEXTA_visualization_4_GMNS.md)

### White Paper and Application

Zhou, X., & Taylor, J. (2014). [DTALite: A queue-based mesoscopic traffic
simulator for fast model evaluation and
calibration.](https://www.tandfonline.com/doi/full/10.1080/23311916.2014.961345)
Cogent Engineering 1.1, 961345.

Marshall, N.L. (2018). [Forecasting the impossible: The status quo of estimating
traffic flows with static traffic assignment and the future of dynamic traffic
assignment.](https://www.sciencedirect.com/science/article/pii/S2210539517301232)
Research in Transportation Business & Management 29, 85-92.

### Youtube Teaching Videos on Use of DTALite/NEXTA Packages

[NeXTA/DTALite Workshop Webinar](https://www.youtube.com/channel/UCUHlqojCQ4f7VvqroUhbaFA) by Jeff
Taylor

### [Mini-Lesson](https://youtu.be/rorZAhNNOf0) on the Algorithmic Details

What is the best way to learn dynamic traffic simulation and network assignment 
for a beginner? Do you want to integrate a powerful traffic simulator in your deep
learning framework? We would like to offer a collaborative learning experience 
through 500 lines of Python codes and real-life data sets. This is part of our 
mini-lessons through teaching dialog.

### Additional References

**1. Parallel computing algorithms**

Qu, Y., & Zhou, X. (2017). Large-scale dynamic transportation network simulation: 
A space-time-event parallel computing approach. Transportation Research Part C: 
Emerging technologies, 75, 1-16.

**2. OD demand estimation**

Lu, C. C., Zhou, X., & Zhang, K. (2013). Dynamic origin–destination demand flow 
estimation under congested traffic conditions.Transportation Research Part C:
Emerging Technologies, 34, 16-37.

**3. Simplified emission estimation model**

Zhou, X., Tanvir, S., Lei, H., Taylor, J., Liu, B., Rouphail, N. M., & Frey, H. C. (2015). 
Integrating a simplified emission estimation model and mesoscopic dynamic traffic 
simulator to efficiently evaluate emission impacts of traffic management strategies.
Transportation Research Part D: Transport and Environment, 37, 123-136.

**4. Eco-system optimal time-dependent flow assignment**

Lu, C. C., Liu, J., Qu, Y., Peeta, S., Rouphail, N. M., & Zhou, X. (2016). Eco-system 
optimal time-dependent flow assignment in a congested network. Transportation Research
Part B: Methodological, 94, 217-239.

**5. Transportation-induced population exposure assessment**

Vallamsundar, S., Lin, J., Konduri, K., Zhou, X., & Pendyala, R. M. (2016). A 
comprehensive modeling framework for transportation-induced population exposure 
assessment. Transportation Research Part D: Transport and Environment, 46, 94-113.

**6. Integrated ABM and DTA**

Xiong, C., Shahabi, M., Zhao, J., Yin, Y., Zhou, X., & Zhang, L. (2020). An 
integrated and personalized traveler information and incentive scheme for energy
efficient mobility systems. Transportation Research Part C: Emerging Technologies, 
113, 57-73.

**7. State-wide transportation modeling**

Zhang, L. (2017). Maryland SHRP2 C10 Implementation Assistance – MITAMS: Maryland 
Integrated Analysis Modeling System, Maryland State Highway Administration.

**8. Workzone applications**

Schroeder, B, et al. (2014). Work zone traffic analysis & impact assessment. FHWA/NC/2012-36. 
North Carolina. Dept. of Transportation. Research and Analysis Group.

> NeXTA Network Visualization
![NeXTA Network Visualization](docs/media/nexta_network.png)

> NeXTA Traffic Visualization
![NeXTA Traffic Visualization](docs/media/nexta_traffic-visualization.png)

> Applications of DTALite/NeXTA in the United States
![Applications of DTALite/NeXTA in the United States](docs/media/dtalite_applications.png)