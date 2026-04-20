package cmd

import (
    "os"
    "github.com/spf13/cobra"
    "github.com/GenomEvo/EvoCG/internal/cg"
    "github.com/GenomEvo/EvoCG/internal/bestpicker"
    "github.com/GenomEvo/EvoCG/internal/batchprot"
    "github.com/GenomEvo/EvoCG/internal/combsingle"
    "github.com/GenomEvo/EvoCG/internal/lenext"
    "github.com/GenomEvo/EvoCG/internal/mutbest"
    "github.com/GenomEvo/EvoCG/internal/protacc"
    "github.com/GenomEvo/EvoCG/internal/simcov"
)

func shiftArgs() {
    // Shift os.Args so that os.Args[1] is the first argument after the subcommand
    if len(os.Args) > 1 {
        newArgs := []string{os.Args[0]}
        newArgs = append(newArgs, os.Args[2:]...)
        os.Args = newArgs
    }
}

var cgCmd = &cobra.Command{
    Use:   "cg",
    Short: "Run the original EvoCG2.0 cg module",
    Run: func(cmd *cobra.Command, args []string) {
        shiftArgs()
        cg.RunCG()
    },
}

var bestpickerCmd = &cobra.Command{
    Use:   "bestpicker",
    Short: "Run bestpicker",
    Run: func(cmd *cobra.Command, args []string) {
        shiftArgs()
        bestpicker.RunBestPicker()
    },
}

var batchProtCmd = &cobra.Command{
    Use:   "batchprot",
    Short: "Run batchProtAcc2locTag",
    Run: func(cmd *cobra.Command, args []string) {
        shiftArgs()
        batchprot.RunBatchProt()
    },
}

var combSingleCmd = &cobra.Command{
    Use:   "combsingle",
    Short: "Run combSingleGeneOrtho",
    Run: func(cmd *cobra.Command, args []string) {
        shiftArgs()
        combsingle.RunCombSingle()
    },
}

var lenExtCmd = &cobra.Command{
    Use:   "lenext",
    Short: "Run lenExt",
    Run: func(cmd *cobra.Command, args []string) {
        shiftArgs()
        lenext.RunLenExt()
    },
}

var mutBestCmd = &cobra.Command{
    Use:   "mutbest",
    Short: "Run mutbest",
    Run: func(cmd *cobra.Command, args []string) {
        shiftArgs()
        mutbest.RunMutBest()
    },
}

var protAccCmd = &cobra.Command{
    Use:   "protacc",
    Short: "Run protAcc2locTag",
    Run: func(cmd *cobra.Command, args []string) {
        shiftArgs()
        protacc.RunProtAcc()
    },
}

var simcovFilterCmd = &cobra.Command{
    Use:   "simcov",
    Short: "Run simcovfilter",
    Run: func(cmd *cobra.Command, args []string) {
        shiftArgs()
        simcov.RunSimCov()
    },
}

func init() {
    rootCmd.AddCommand(cgCmd)
    rootCmd.AddCommand(bestpickerCmd)
    rootCmd.AddCommand(batchProtCmd)
    rootCmd.AddCommand(combSingleCmd)
    rootCmd.AddCommand(lenExtCmd)
    rootCmd.AddCommand(mutBestCmd)
    rootCmd.AddCommand(protAccCmd)
    rootCmd.AddCommand(simcovFilterCmd)
}
