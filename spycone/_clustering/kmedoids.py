# this script is derived from sklearn-extra packages KMedoids

import numpy as np
import warnings

from sklearn.base import BaseEstimator, ClusterMixin, TransformerMixin

from . import similarities as sim

class KMedoids(BaseEstimator, ClusterMixin, TransformerMixin):
    def __init__(
        self,
        n_clusters=8,
        metric="precomputed",
        max_iter=300,
        random_state=None,
        method="pam"
    ):
        self.n_clusters = n_clusters
        self.metric = metric
        self.max_iter = max_iter
        self.random_state = random_state
        self.method = method

    def _initialize_medoids(self, D, n_clusters):
        """define initial medoids:
        heuristic""" 

        medoids = np.random.choice(D.shape[0], n_clusters)
        return medoids


    def _compute_inertia(self, distances):
        """Compute inertia of new samples. Inertia is defined as the sum of the
        sample distances to closest cluster centers.
        Parameters
        ----------
        distances : {array-like, sparse matrix}, shape=(n_samples, n_clusters)
            Distances to cluster centers.
        Returns
        -------
        Sum of sample distances to closest cluster centers.
        """

        # Define inertia as the sum of the sample-distances
        # to closest cluster centers
        inertia = np.sum(np.min(distances, axis=1))

        return inertia

    def transform(self, X):
        """Transforms X to cluster-distance space.
        Parameters
        ----------clu = clustering(allcovdata, algorithm='kmedoids', metric='correlation',BioNetwork=bionet, composite=False)
c = clu.find_clusters()

%matplotlib inline
cluster_heatmap(clu, standard_scale=0)
            X transformed in the new space of distances to cluster centers.
        """
    
        if self.metric == "precomputed":
            return X[:, self.medoid_indices_]

    def _update_medoid_idxs_in_place(self, D, labels, medoid_idxs):
        """In-place update of the medoid indices"""

        # Update the medoids for each cluster
        for k in range(self.n_clusters):
            # Extract the distance matrix between the data points
            # inside the cluster k
            cluster_k_idxs = np.where(labels == k)[0]

            if len(cluster_k_idxs) == 0:
                warnings.warn(
                    "Cluster {k} is empty! "
                    "self.labels_[self.medoid_indices_[{k}]] "
                    "may not be labeled with "
                    "its corresponding cluster ({k}).".format(k=k)
                )
                continue

            in_cluster_distances = D[
                cluster_k_idxs, cluster_k_idxs[:, np.newaxis]
            ]

            # Calculate all costs from each point to all others in the cluster
            in_cluster_all_costs = np.sum(in_cluster_distances, axis=1)

            min_cost_idx = np.argmin(in_cluster_all_costs)
            min_cost = in_cluster_all_costs[min_cost_idx]
            curr_cost = in_cluster_all_costs[
                np.argmax(cluster_k_idxs == medoid_idxs[k])
            ]

            # Adopt a new medoid if its distance is smaller then the current
            if min_cost < curr_cost:
                medoid_idxs[k] = cluster_k_idxs[min_cost_idx]

    def _compute_cost(self, D, medoid_idxs):
        """ Compute the cose for a given configuration of the medoids"""
        return self._compute_inertia(D[:, medoid_idxs])

    def fit(self, X):
        """
        x: distance array
        """
        if self.n_clusters > X.shape[0]:
            raise ValueError(
                "The number of medoids (%d) must be less "
                "than the number of samples %d."
                % (self.n_clusters, X.shape[0])
            )

        medoid_idxs = self._initialize_medoids(X, self.n_clusters)

        #compute the distance to the first and second closest points among medoids
        Djs, Ejs = np.sort(X[medoid_idxs], axis=0)[[0,1]]


        for self.n_iter_ in range(0, self.max_iter):
            old_medoid_idxs = np.copy(medoid_idxs)
            labels = np.argmin(X[medoid_idxs,:], axis=0)

            if self.method == "alternate":
                # Update medoids with the new cluster indices
                self._update_medoid_idxs_in_place(D, labels, medoid_idxs)
            
            elif self.method == "pam":
                not_medoid_idxs = np.delete(np.arange(len(X)), medoid_idxs)
                optimal_swap = sim._compute_optimal_swap(X,
                                                        medoid_idxs.astype(np.intc),
                                                        not_medoid_idxs.astype(np.intc),
                                                        Djs,
                                                        Ejs,
                                                        self.n_clusters)

                if optimal_swap is not None:
                    i,j, _ = optimal_swap
                    medoid_idxs[medoid_idxs == i] = j

                    #update
                    Djs, Ejs = np.sort(X[medoid_idxs], axis=0)[[0,1]]

            if np.all(old_medoid_idxs == medoid_idxs):
                break
            elif self.n_iter_ == self.max_iter - 1:
                warnings.warn(
                    "Maximum number of iteration reached before "
                    "convergence. Consider increasing max_iter to "
                    "improve the fit."
                )

            # Expose labels_ which are the assignments of
            # the training data to clusters
            self.labels_ = np.argmin(X[medoid_idxs, :], axis=0)
            self.medoid_indices_ = medoid_idxs
            self.inertia_ = self._compute_inertia(self.transform(X))

            # Return self to enable method chaining
            return self