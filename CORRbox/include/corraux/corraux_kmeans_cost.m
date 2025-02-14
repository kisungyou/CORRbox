function cost = corraux_kmeans_cost(variation, label)
    
    % when a label vector and a variation vector is given, compute the cost
    % objective value for the standard k-means algorithm. 

    ulabel = sort(unique(label));
    K      = length(ulabel);
    
    cost   = 0;
    for k=1:K
        cost = cost + variation(k)*sum(label==ulabel(k));
    end
end