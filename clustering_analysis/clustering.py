
def find_peptide_clusters(matrix):
    clusters = []
    for i in range(len(matrix)):
        current_cluster = None
        for cluster in clusters:
            if i in cluster:
                current_cluster = cluster
        if current_cluster is None:
            current_cluster = set()
            clusters.append(current_cluster)
        current_cluster.add(i)
        for j in range(len(matrix[i])):
            if matrix[i][j] == 1:
                current_cluster.add(j)

    clusters = merge_clusters(clusters)

    assert clusters is not None
    assert all(isinstance(c, set) for c in clusters)

    total_clustered = sum(len(c) for c in clusters)
    unique_clustered = len(set().union(*clusters))
    assert total_clustered == unique_clustered  # no peptide in two clusters
    assert unique_clustered == 30               # every peptide assigned

    return clusters


def merge_clusters(clusters):
    changed = True

    while changed:
        changed = False
        i = 0

        while i < len(clusters):
            j = i + 1

            while j < len(clusters):
                if clusters[i] & clusters[j]:
                    clusters[i] |= clusters[j]
                    del clusters[j]
                    changed = True
                else:
                    j += 1

            i += 1

    return clusters
