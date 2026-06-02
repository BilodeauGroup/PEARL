
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

    i = 0
    while i < len(clusters):
        j = i + 1
        while j < len(clusters):
            if clusters[i] & clusters[j]:
                clusters[i] |= clusters[j]
                del clusters[j]
            else:
                j += 1
        i += 1
    return clusters
