classdef RestrictClusters < handle
% RestrictClusters Class for
%
%  RestrictClusters Properties:
%       isi_viol        - ISI violation threshold (default = 0.15)
%       isolation_dist  - isolation distance threshold (default = 20)
%       amp_cutoff      - amplitude cutoff threshold (default = 0.1)
%       max_drift       - max. drift allowed (default = 80)
%       ctl             - instance of RC2Preprocess
%       probe_id        - string with probe recording ID
%
%  RestrictClusters Methods:
%       curation_table          - automated Bombcell + metric curation as an
%                                 editable per-cluster table (with `keep` column)
%       curation_mua_table      - same, for MUA clusters
%       restrict_metrics_table  - (legacy) metrics.csv restricted to clusters
%                                 passing quality criteria AND Bombcell 'good'
%       restrict_mua_metrics_table - (legacy) MUA variant

    properties
        
        isi_viol        = 0.15
        isolation_dist  = 20
        amp_cutoff      = 0.1
        max_drift       = 80
    end
    
    properties (SetAccess = private)
        
        ctl
        probe_id
    end
    
    
    methods
        
        function obj = RestrictClusters(ctl, probe_id)
        %%RestrictClusters
        %
        %   RestrictClusters(CTL, PROBE_ID). CTL is instance of
        %   RC2Preprocess and PROBE_ID is a string with the probe recording
        %   ID.
        
            obj.ctl = ctl;
            obj.probe_id = probe_id;
        end
        
        
        
        function new_tbl = restrict_metrics_table(obj)
        %%restrict_metrics_table Loads the metrics.csv and restricts it to
        % clusters which satisfy certain quality criteria
        %
        %   TABLE = restrict_metrics_table() restricts the metrics.csv to
        %   clusters which satisfy quality of criteria. 
        %       < `isi_viol`
        %       > `isolation_distance`
        %       < `amp_cutoff`
        %       < `max_drift`
        
            metrics = obj.ctl.load.metrics_csv(obj.probe_id);
            cluster_groups = obj.ctl.load.cluster_groups(obj.probe_id);
            
            % NaN-tolerant comparisons: if a metric could not be computed
            % (NaN), the cluster is not excluded on that criterion alone.
            % This is necessary for KS4 where isolation_distance is NaN
            % for most clusters (quality_metrics module limitation with
            % KS4 PC feature format).
            iso_ok   = isnan(metrics.isolation_distance) | metrics.isolation_distance > obj.isolation_dist;
            drift_ok = isnan(metrics.max_drift)          | metrics.max_drift          < obj.max_drift;

            idx = metrics.isi_viol < obj.isi_viol & ...
                  iso_ok & ...
                  metrics.amplitude_cutoff < obj.amp_cutoff & ...
                  drift_ok;
            
            good_idx = false(size(metrics, 1), 1);
              
            for ii = 1 : size(metrics, 1)
                
                clust_idx = find(metrics.cluster_id(ii) == cluster_groups.cluster_id);
                
                if isempty(clust_idx)
                    continue
                end
                
                good_idx(ii) = strcmp(cluster_groups.group(clust_idx), 'good');
            end
            
            idx = idx & good_idx;
            
            new_tbl = metrics(idx, :);
        end
        
        
        
        function new_tbl = restrict_mua_metrics_table(obj)
        %%restrict_mua_metrics_table Loads the metrics.csv and restricts it
        %%to clusters which satisfy certain quality criteria for our MUA
        %%clusters
        %
        %   TABLE = restrict_mua_metrics_table() restricts the metrics.csv
        %   to clusters which satisfy quality of criteria. 
        %       < `amp_cutoff`
        %       < `max_drift`
        
            metrics = obj.ctl.load.metrics_csv(obj.probe_id);
            
            % remove heading with Var1
            if ismember('Var1', metrics.Properties.VariableNames)
                metrics.Var1 = [];
            end
            
            drift_ok = isnan(metrics.max_drift) | metrics.max_drift < obj.max_drift;
            idx = metrics.amplitude_cutoff < obj.amp_cutoff & drift_ok;

            new_tbl = metrics(idx, :);
        end



        function tbl = curation_table(obj)
        %%curation_table Automated Bombcell + quality-metric curation, as an editable table
        %
        %   TABLE = curation_table() returns ONE ROW PER CLUSTER with columns:
        %       cluster_id      - Kilosort cluster ID
        %       bombcell_group  - raw Bombcell label (good/mua/noise/
        %                         non_soma_good/non_soma_mua)
        %       is_non_somatic  - 1 if Bombcell flagged the unit as axonal/
        %                         dendritic (bombcell_group starts with
        %                         'non_soma_'), 0 otherwise
        %       passes_metrics  - 1 if the cluster passes the quality thresholds
        %                         (isi_viol, isolation_distance, amp_cutoff, max_drift)
        %       keep_pipeline   - the automated decision, 1 = keep: Bombcell
        %                         'good' (somatic only) AND passes_metrics.
        %                         Left UNTOUCHED by hand -- a fixed record of
        %                         what the pipeline alone decided, so a later
        %                         hand-edit of `keep` can always be compared
        %                         back against it.
        %       keep            - starts identical to keep_pipeline. Non-somatic
        %                         units default to 0 even if labelled
        %                         'non_soma_good' -- inspect is_non_somatic/
        %                         bombcell_group and hand-edit keep=1 to
        %                         include them (e.g. for axonal signal
        %                         analyses).
        %
        %   Only the `keep` column is consumed downstream
        %   (create_selected_clusters_txt), so a user may hand-edit it to force
        %   any cluster in (keep=1) or out (keep=0) before formatting -- the
        %   selection is otherwise fully automated. Do not edit keep_pipeline;
        %   it is regenerated from scratch every time create_check_clusters_csv
        %   runs and exists only so a hand-edited keep can be checked against
        %   the original automated decision later.

            metrics        = obj.ctl.load.metrics_csv(obj.probe_id);
            cluster_groups = obj.ctl.load.cluster_groups(obj.probe_id);

            % NaN-tolerant quality-threshold pass (same rules as restrict_metrics_table)
            iso_ok   = isnan(metrics.isolation_distance) | metrics.isolation_distance > obj.isolation_dist;
            drift_ok = isnan(metrics.max_drift)          | metrics.max_drift          < obj.max_drift;
            passes_metrics = metrics.isi_viol < obj.isi_viol & ...
                             iso_ok & ...
                             metrics.amplitude_cutoff < obj.amp_cutoff & ...
                             drift_ok;

            [bombcell_group, is_good, is_non_somatic] = ...
                obj.lookup_bombcell_group(metrics.cluster_id, cluster_groups, 'good');

            keep_pipeline = double(is_good & ~is_non_somatic & passes_metrics);
            keep          = keep_pipeline;

            tbl = table(metrics.cluster_id, bombcell_group, double(is_non_somatic), ...
                double(passes_metrics), keep_pipeline, keep, ...
                'VariableNames', {'cluster_id', 'bombcell_group', 'is_non_somatic', 'passes_metrics', 'keep_pipeline', 'keep'});
        end



        function tbl = curation_mua_table(obj)
        %%curation_mua_table Automated MUA curation, as an editable table
        %
        %   TABLE = curation_mua_table() returns ONE ROW PER CLUSTER with the
        %   same columns as curation_table (including keep_pipeline, the fixed
        %   record of the automated decision -- see curation_table), but
        %   `keep`/`keep_pipeline` default to 1 when the Bombcell label is
        %   'mua' (somatic only) AND the cluster passes the (looser) MUA
        %   metric thresholds (amp_cutoff, max_drift). Non-somatic units
        %   default to keep=0 -- see curation_table for how to override.
        %   Hand-edit `keep` (never keep_pipeline) to override before
        %   create_selected_mua_clusters_txt.

            metrics        = obj.ctl.load.metrics_csv(obj.probe_id);
            cluster_groups = obj.ctl.load.cluster_groups(obj.probe_id);

            if ismember('Var1', metrics.Properties.VariableNames)
                metrics.Var1 = [];
            end

            drift_ok = isnan(metrics.max_drift) | metrics.max_drift < obj.max_drift;
            passes_metrics = metrics.amplitude_cutoff < obj.amp_cutoff & drift_ok;

            [bombcell_group, is_mua, is_non_somatic] = ...
                obj.lookup_bombcell_group(metrics.cluster_id, cluster_groups, 'mua');

            keep_pipeline = double(is_mua & ~is_non_somatic & passes_metrics);
            keep          = keep_pipeline;

            tbl = table(metrics.cluster_id, bombcell_group, double(is_non_somatic), ...
                double(passes_metrics), keep_pipeline, keep, ...
                'VariableNames', {'cluster_id', 'bombcell_group', 'is_non_somatic', 'passes_metrics', 'keep_pipeline', 'keep'});
        end
    end



    methods (Static = true)

        function [group, is_label, is_non_somatic] = lookup_bombcell_group(cluster_ids, cluster_groups, label)
        %%lookup_bombcell_group Map each cluster_id to its Bombcell group string
        %
        %   [GROUP, IS_LABEL, IS_NON_SOMATIC] = lookup_bombcell_group(CLUSTER_IDS, CLUSTER_GROUPS, LABEL)
        %   returns GROUP (cell array of the raw Bombcell label per cluster, ''
        %   if not found), IS_LABEL (logical, true where the underlying somatic
        %   label -- 'good' or 'mua' -- equals LABEL, whether or not the unit
        %   is non-somatic), and IS_NON_SOMATIC (logical, true for
        %   'non_soma_good' / 'non_soma_mua', i.e. axonal/dendritic units;
        %   see bombcell_label_units split_non_somatic_good_mua).

            n = numel(cluster_ids);
            group          = repmat({''}, n, 1);
            is_label       = false(n, 1);
            is_non_somatic = false(n, 1);
            for ii = 1 : n
                idx = find(cluster_ids(ii) == cluster_groups.cluster_id, 1);
                if isempty(idx)
                    continue
                end
                group{ii} = cluster_groups.group{idx};
                is_non_somatic(ii) = startsWith(group{ii}, 'non_soma_');
                if is_non_somatic(ii)
                    underlying_label = erase(group{ii}, 'non_soma_');
                else
                    underlying_label = group{ii};
                end
                is_label(ii) = strcmp(underlying_label, label);
            end
        end
    end
end
