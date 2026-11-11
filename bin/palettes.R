# palettes.R — shared colour helpers for bulkAnnex pipeline scripts.
#
# ba_distinct_pal() returns a *deterministic* set of well-separated colours for a
# set of categorical levels, using only base R (grDevices) — no new packages.
# Levels are sorted internally so that the same set of levels always maps to the
# same colours, independent of call order. This keeps the pipeline PNGs and the
# dashboard in sync (a given sample/condition always gets the same colour).
#
# Strategy: a fixed, ordered panel (`.ba_panel`) is indexed from the front, so the
# most-distinguishable colours are used first and redundant colours are only
# reached when there are very many categories:
#   1. BuenColors "corona"  (30)  — tab10 + tab20 pastels, the preferred look
#   2. RColorBrewer Dark2   (8)
#   3. RColorBrewer Set1    (9)
#   4. RColorBrewer Set2    (8)   — last resort (pastel, overlaps corona)
#   5. RColorBrewer Set3    (12)  — last resort (pastel, overlaps corona)
# Beyond 67 categories a golden-angle hue wheel is used (no duplicate hex, but no
# palette can make that many categories truly distinguishable).
.ba_panel <- c(
    # corona — 30
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf", "#ad494a", "#e7ba52", "#8ca252", "#756bb1", "#636363", "#aec7e8",
    "#ff9896", "#98df8a", "#ffbb78", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d",
    "#9edae5", "#e7969c", "#e7cb94", "#c7e9c0", "#de9ed6", "#d9d9d9",
    # Dark2 — 8
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
    # Set1 — 9
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
    # Set2 — 8 (last resort)
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
    # Set3 — 12 (last resort)
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
    "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)

ba_distinct_pal <- function(levels) {
    levels <- as.character(levels)
    levels <- sort(unique(levels[!is.na(levels) & nzchar(levels)]))
    n <- length(levels)
    cols <- if (n == 0) {
        character(0)
    } else if (n <= length(.ba_panel)) {
        .ba_panel[seq_len(n)]
    } else {
        grDevices::hcl(h = (seq_len(n) * 137.508) %% 360, c = 55, l = 65)
    }
    stats::setNames(cols, levels)
}
