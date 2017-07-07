[
    {                           // main CLI
        type: "wire-cell",
        data: {
            // fixme: need gen for AnodePlane, best to move that to a 3rd lib
            plugins: ["WireCellGen", "WireCellSigProc", "WireCellSio"],
            apps: ["Omnibus"]
        }
    },
]
