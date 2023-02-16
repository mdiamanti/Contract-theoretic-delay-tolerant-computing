# Contract-theoretic-delay-tolerant-computing
Code for generating An Incentivization Mechanism for Green Computing Continuum of Delay-Tolerant Tasks paper's simulation environment. https://ieeexplore.ieee.org/document/9838752

# Abstract
Capitalizing on the different available computing options across the network, the concept of computing continuum has recently emerged to efficiently manage the exaggerated computation demands of the numerous Internet-of-Things (IoT) users and applications. Nevertheless, the edge computing’s attractiveness to the users, in terms of its reduced incurred time and energy overhead, acts as an impediment in the realization of the envisioned computing continuum. In this paper, recognizing the potential of forwarding delay-tolerant tasks to upper computing layers, we design an incentivization-based mechanism for the offloading users, aiming to shift their preference from the edge to the upper fog computing layer. The corresponding mechanism comprises two stages, in which different models of Contract Theory are adopted. In the first stage, a users-to-edge server contract is formulated to determine the optimal amount of each user’s initially offloaded task at the edge that is allowed to be
further forwarded to the fog, based on the user’s delay tolerance. Subsequently, an edge-to-fog server contract is formulated to account for the edge server’s tradeoff between the local execution and transmission overheads, deriving the most beneficial amount of the users’ tasks that ultimately reaches the fog. The overall mechanism is evaluated via modeling and simulation regarding its operation and efficiency under different scenarios.

# Run simulations
To generate numerical results from the user to the edge, execute results_user_to_edge.m
To generate numerical results from the edge to fog, execute results_edge_to_fog.m
