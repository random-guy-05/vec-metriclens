TASK_METRICS = {
    "T1": {
        "de_score": ("higher", 0.25, "DE gene recovery"),
        "de_direction": ("higher", 0.25, "Change direction"),
        "mmd_u": ("lower", 0.30, "Cell-state distribution"),
        "variogram": ("lower", 0.20, "Gene-gene co-variation"),
    },
    "T2": {
        "de_score": ("higher", 0.125, "Expression change"),
        "de_direction": ("higher", 0.125, "Expression change"),
        "mmd_u": ("lower", 0.15, "Cell-state distribution"),
        "variogram": ("lower", 0.10, "Cell-state distribution"),
        "d2_shape": ("lower", 1 / 12, "Tissue shape and growth scale"),
        "occupancy_dice": ("higher", 1 / 12, "Tissue shape and growth scale"),
        "scale_log_ratio": ("zero", 1 / 12, "Tissue shape and growth scale"),
        "neighborhood_mmd": ("lower", 0.25, "Local spatial organisation"),
    },
    "T3": {
        "de_score": ("higher", 0.30, "Response gene recovery"),
        "de_direction": ("higher", 0.25, "Response direction"),
        "severity_slope": ("zero", 0.25, "Response magnitude"),
        "mmd_u": ("lower", 0.12, "Cell-state distribution"),
        "variogram": ("lower", 0.08, "Cell-state distribution"),
    },
}

ANCHORS = {
    "T1:val": {
        "task": "T1",
        "metrics": {
            "de_score": (0.0, 0.8464),
            "de_direction": (0.0, 0.7901),
            "mmd_u": (0.08359, 0.00406),
            "variogram": (0.005219, 0.000158),
        },
    },
    "T2:embryo:val_interp": {
        "task": "T2",
        "metrics": {
            "de_score": (0.0, 0.8413),
            "de_direction": (0.0, 0.9185),
            "mmd_u": (0.08571, 0.00307),
            "variogram": (0.053533, 0.001264),
            "d2_shape": (0.05306, 0.00268),
            "occupancy_dice": (0.7047, 0.7661),
            "scale_log_ratio": (-0.3063, 0.0053),
            "neighborhood_mmd": (0.21105, 0.01128),
        },
    },
    "T2:heart:val_extrap": {
        "task": "T2",
        "metrics": {
            "de_score": (0.0, 0.9420),
            "de_direction": (0.0, 0.9915),
            "mmd_u": (0.02455, 0.00011),
            "variogram": (0.031712, 0.000696),
            "d2_shape": (0.00933, 0.00342),
            "occupancy_dice": (0.7747, 0.9465),
            "scale_log_ratio": (-0.2680, 0.0074),
            "neighborhood_mmd": (0.07184, 0.0004),
        },
    },
    "T2:heart:val_interp": {
        "task": "T2",
        "metrics": {
            "de_score": (0.0, 0.8182),
            "de_direction": (0.0, 0.9553),
            "mmd_u": (0.0207, 0.00087),
            "variogram": (0.023828, 0.001021),
            "d2_shape": (0.03079, 0.00412),
            "occupancy_dice": (0.6748, 0.8796),
            "scale_log_ratio": (0.3284, -0.0119),
            "neighborhood_mmd": (0.05728, 0.00432),
        },
    },
    "T3:gata4": {
        "task": "T3",
        "metrics": {
            "de_score": (0.0, 0.8696),
            "de_direction": (0.0, 0.9247),
            "severity_slope": (-6.9078, -0.0088),
            "mmd_u": (0.03141, 0.00039),
            "variogram": (0.032212, 0.000737),
        },
    },
}
