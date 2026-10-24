import * as THREE from 'three';

export const SIGNAL_COLOR = '#43c6d5';
export const MAX_SIGNAL_BLEND = .72;

// A display-only color wash over the original shaded material. Nothing is cut,
// made transparent, deformed, or turned into a predicted spatial drug field.
export function attachAmountSignal(root) {
  return (root.userData.signalMaterials || []).map(material => {
    const amount = { value: 0 };
    const signalColor = { value: new THREE.Color(SIGNAL_COLOR) };
    material.onBeforeCompile = shader => {
      shader.uniforms.cabombAmount = amount;
      shader.uniforms.cabombSignal = signalColor;
      shader.fragmentShader = `uniform float cabombAmount;\nuniform vec3 cabombSignal;\n${shader.fragmentShader}`;
      shader.fragmentShader = shader.fragmentShader.replace('#include <opaque_fragment>', `
        // Keep edge shading and structural texture legible beneath the wash.
        float cabombFacing = abs(dot(normalize(normal), normalize(vViewPosition)));
        float cabombMask = mix(0.30, 1.0, smoothstep(0.16, 0.78, cabombFacing));
        float cabombLuma = dot(outgoingLight, vec3(0.2126, 0.7152, 0.0722));
        vec3 cabombWash = cabombSignal * (0.23 + 0.95 * cabombLuma);
        outgoingLight = mix(outgoingLight, cabombWash, cabombAmount * ${MAX_SIGNAL_BLEND.toFixed(2)} * cabombMask);
        #include <opaque_fragment>
      `);
    };
    material.customProgramCacheKey = () => 'cabomb-original-shape-wash-v1';
    material.needsUpdate = true;
    return { amount, material };
  });
}
