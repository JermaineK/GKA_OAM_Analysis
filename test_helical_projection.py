import numpy as np

from spiral_decoder_experiment import helical_projection, theta_image_coords


def test_helical_projection_matches_signed_azimuthal_phase() -> None:
    """
    Synthetic field with pure azimuthal phase exp(i*ell*theta) should project back to the same ell.
    This guards against handedness flips from image-coordinate conventions.
    """
    h = w = 128
    center = (w / 2.0, h / 2.0)
    yy, xx = np.mgrid[0:h, 0:w]
    theta = theta_image_coords(yy, xx, center[1], center[0])
    ells = (-3, -2, -1, 1, 2, 3)

    for ell in ells:
        field = np.exp(1j * ell * theta).astype(np.complex64)
        proj = helical_projection(
            field,
            center=center,
            ells=ells,
            r_inner_frac=0.08,
            r_outer_frac=0.4,
        )
        best_ell, best_val = max(
            (
                (int(k.split("proj_ell_")[-1]), v)
                for k, v in proj.items()
                if k.startswith("proj_ell_")
            ),
            key=lambda kv: kv[1],
        )
        assert best_ell == ell, f"Expected ell={ell}, got {best_ell} (|c|={best_val:.3f})"
        assert best_val > 0.8, f"Projection magnitude too small for ell={ell}: {best_val:.3f}"
