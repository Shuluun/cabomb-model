# Result preview

## Delivery-model structure

![Circulation and targeting compartments](../assets/circulation_targeting_diagram.svg)

The amount-based model separates circulation and organ passage from tumor
entry, binding and internalization. OCN and YSA act after tumor-region arrival.

## Synthetic forecasting example

![Early-only forecasts on synthetic data](../assets/synthetic_forecast.png)

The vertical line marks the end of the observed early window. Forecasts use
those observations; later values are scored separately. This noise-free example
illustrates the runnable workflow. Experimental comparisons are discussed in
the research notes.

## Release-waiting comparison

![Assumed waiting-time responses](../assets/bridge_response.png)

Each panel uses its own y-axis scale. Three assumed waiting shapes share a
four-hour mean and the same model-generated internalization input for that
design. Their early timing differs, while late responses approach one another.
Units are normalized carrier equivalents. The input and output CSVs are saved
separately when the commands run.
