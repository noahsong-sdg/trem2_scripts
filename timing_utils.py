#!/usr/bin/env python3
"""
Timing utilities for tracking performance metrics in the UniDock HTS workflow.
Provides consistent timing, progress tracking, and performance reporting across all scripts.
"""

import time
import json
import os
from datetime import datetime, timedelta
from pathlib import Path

class TimingTracker:
    """
    A comprehensive timing tracker for HTS workflows.
    Tracks overall runtime, step-by-step timing, and calculates performance metrics.
    """
    
    def __init__(self, script_name, log_dir=None):
        self.script_name = script_name
        self.start_time = time.time()
        self.step_times = {}
        self.step_start_time = None
        self.current_step = None
        self.ligands_processed = 0
        self.total_ligands = 0
        
        # Set up logging directory
        if log_dir is None:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            log_dir = os.path.join(script_dir, "../results/timing_logs")
        
        self.log_dir = log_dir
        os.makedirs(self.log_dir, exist_ok=True)
        
        # Create timestamped log file
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_file = os.path.join(self.log_dir, f"{script_name}_{timestamp}.json")
        
        print(f"🕐 Starting {script_name} at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"📊 Timing data will be saved to: {self.log_file}")
    
    def start_step(self, step_name, total_items=None):
        """Start timing a specific step."""
        if self.current_step:
            self.end_step()
        
        self.current_step = step_name
        self.step_start_time = time.time()
        
        if total_items:
            self.total_ligands = total_items
            print(f"\n🚀 Starting step: {step_name} ({total_items} items to process)")
        else:
            print(f"\n🚀 Starting step: {step_name}")
    
    def end_step(self):
        """End timing the current step."""
        if not self.current_step:
            return
        
        elapsed = time.time() - self.step_start_time
        self.step_times[self.current_step] = elapsed
        
        # Calculate performance metrics for this step
        if self.ligands_processed > 0:
            ligands_per_minute = (self.ligands_processed / elapsed) * 60
            avg_time_per_ligand = elapsed / self.ligands_processed
            
            print(f"✅ Completed step: {self.current_step}")
            print(f"   Duration: {self.format_duration(elapsed)}")
            print(f"   Processed: {self.ligands_processed:,} ligands")
            print(f"   Rate: {ligands_per_minute:.1f} ligands/minute")
            print(f"   Average: {avg_time_per_ligand:.3f} seconds/ligand")
        else:
            print(f"✅ Completed step: {self.current_step} in {self.format_duration(elapsed)}")
        
        self.current_step = None
        self.step_start_time = None
        self.ligands_processed = 0
    
    def update_progress(self, processed_count, step_increment=100):
        """Update progress within a step."""
        self.ligands_processed = processed_count
        
        if processed_count > 0 and processed_count % step_increment == 0:
            if self.step_start_time:
                elapsed = time.time() - self.step_start_time
                rate = (processed_count / elapsed) * 60
                
                # Estimate time remaining
                if self.total_ligands > 0:
                    remaining = self.total_ligands - processed_count
                    eta_seconds = remaining / (processed_count / elapsed) if processed_count > 0 else 0
                    eta = self.format_duration(eta_seconds)
                    progress_pct = (processed_count / self.total_ligands) * 100
                    print(f"   Progress: {processed_count:,}/{self.total_ligands:,} ({progress_pct:.1f}%) | "
                          f"Rate: {rate:.1f}/min | ETA: {eta}")
                else:
                    print(f"   Progress: {processed_count:,} processed | Rate: {rate:.1f}/min")
    
    def finish(self):
        """Finish timing and generate final report."""
        if self.current_step:
            self.end_step()
        
        total_elapsed = time.time() - self.start_time
        
        # Generate comprehensive report
        report = {
            "script_name": self.script_name,
            "start_time": datetime.fromtimestamp(self.start_time).isoformat(),
            "end_time": datetime.now().isoformat(),
            "total_duration_seconds": total_elapsed,
            "total_duration_formatted": self.format_duration(total_elapsed),
            "step_timings": {},
            "performance_metrics": {}
        }
        
        # Add step timing details
        for step, duration in self.step_times.items():
            report["step_timings"][step] = {
                "duration_seconds": duration,
                "duration_formatted": self.format_duration(duration),
                "percentage_of_total": (duration / total_elapsed) * 100
            }
        
        # Calculate overall performance metrics
        total_ligands_all_steps = sum(getattr(self, f'ligands_step_{i}', 0) for i in range(10))
        if hasattr(self, 'final_ligand_count') and self.final_ligand_count > 0:
            ligands_per_minute = (self.final_ligand_count / total_elapsed) * 60
            avg_time_per_ligand = total_elapsed / self.final_ligand_count
            
            report["performance_metrics"] = {
                "total_ligands_processed": self.final_ligand_count,
                "ligands_per_minute": ligands_per_minute,
                "average_seconds_per_ligand": avg_time_per_ligand,
                "estimated_time_for_1M_ligands": self.format_duration((1_000_000 / ligands_per_minute) * 60),
                "estimated_time_for_10M_ligands": self.format_duration((10_000_000 / ligands_per_minute) * 60)
            }
        
        # Save detailed report
        with open(self.log_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Print summary
        print(f"\n🏁 {self.script_name} completed!")
        print(f"📊 Total runtime: {self.format_duration(total_elapsed)}")
        
        if "performance_metrics" in report and report["performance_metrics"]:
            metrics = report["performance_metrics"]
            print(f"🚀 Performance Summary:")
            print(f"   Ligands processed: {metrics.get('total_ligands_processed', 0):,}")
            print(f"   Processing rate: {metrics.get('ligands_per_minute', 0):.1f} ligands/minute")
            print(f"   Average time per ligand: {metrics.get('average_seconds_per_ligand', 0):.3f} seconds")
            print(f"\n📈 Scale Estimates:")
            print(f"   1M ligands would take: {metrics.get('estimated_time_for_1M_ligands', 'N/A')}")
            print(f"   10M ligands would take: {metrics.get('estimated_time_for_10M_ligands', 'N/A')}")
        else:
            print(f"🚀 Performance Summary: No successful ligand processing completed")
        
        print(f"📋 Detailed timing report saved to: {self.log_file}")
        
        return report
    
    def set_final_ligand_count(self, count):
        """Set the final count of successfully processed ligands."""
        self.final_ligand_count = count
    
    @staticmethod
    def format_duration(seconds):
        """Format duration in a human-readable way."""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            minutes = int(seconds // 60)
            secs = seconds % 60
            return f"{minutes}m {secs:.1f}s"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            secs = seconds % 60
            return f"{hours}h {minutes}m {secs:.1f}s"

def load_timing_reports(log_dir=None):
    """Load and analyze all timing reports."""
    if log_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        log_dir = os.path.join(script_dir, "../results/timing_logs")
    
    if not os.path.exists(log_dir):
        print(f"No timing logs found in {log_dir}")
        return []
    
    reports = []
    for file in Path(log_dir).glob("*.json"):
        try:
            with open(file, 'r') as f:
                report = json.load(f)
                reports.append(report)
        except Exception as e:
            print(f"Error loading {file}: {e}")
    
    return sorted(reports, key=lambda x: x.get('start_time', ''), reverse=True)

def print_performance_summary():
    """Print a summary of recent performance metrics."""
    reports = load_timing_reports()
    
    if not reports:
        print("No timing reports found.")
        return
    
    print("\n📊 Recent Performance Summary:")
    print("=" * 80)
    
    for report in reports[:5]:  # Show last 5 runs
        script = report['script_name']
        duration = report['total_duration_formatted']
        start_time = report.get('start_time', 'Unknown').split('T')[0]
        
        print(f"\n{script} ({start_time}):")
        print(f"  Duration: {duration}")
        
        if 'performance_metrics' in report:
            metrics = report['performance_metrics']
            print(f"  Ligands: {metrics['total_ligands_processed']:,}")
            print(f"  Rate: {metrics['ligands_per_minute']:.1f} ligands/min")

if __name__ == "__main__":
    print_performance_summary()
